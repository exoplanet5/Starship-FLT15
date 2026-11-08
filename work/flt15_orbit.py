"""Starship Flight 15 (NET 13 Nov 2026, window 22:45-00:00 UTC): orbit model, nominal ascent/coast/orbit, Super Heavy return.

Plane: Flight 14's fit to the September navigational warnings (ift14_fit: i = 30.58 deg, 275 km circular with J2 secular rates,
node 2.1 deg east of the pad, the ascent steering into the plane). SpaceX: Flight 15 flies the same altitude and inclination.
Timeline: the official Flight 15 flight timeline (../flight-timeline.txt). The Earth-fixed track is the same for every T0.

Ascent: the Flight 14 ground path (cubic Hermite from the pad at rest to the SECO state) re-timed along its arc length with a
shaped distance profile through typical Starship telemetry anchors (Max Q, MECO, SECO), and a height profile that ends on the
coast ellipse in height and climb rate. SECO height (150 km) and downrange (15 deg of arc) are assumptions carried from Flight 14.
Super Heavy: 2-D point mass from the MECO state - coast, 36 s boostback, ballistic fall with US76 drag, 25 s landing burn to the
tower - with the boostback acceleration and pitch solved so the fall reaches the landing-burn start over the tower at T+06:10."""
import sys, numpy as np
from pathlib import Path
from scipy.optimize import brentq, least_squares
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
sys.path.insert(0, '/Users/mickey/sda/starship-ITF-14/work')
import ift14_fit as F                      # Flight 14 plane fit, zones, Kepler helpers
from ift14_atmo import atmo

ROOT = Path('/Users/mickey/launch/Starship-FLT15')
MU, RE, J2, OMEGA_E, EF2 = F.MU, F.RE, F.J2, F.OMEGA_E, F.EF2       # OMEGA_E in rad/min
STARBASE = F.STARBASE                                               # pad used by the Flight 14 fit (lon, lat)
TOWER = (-97.1545, 25.9965)                                         # Pad 2 catch tower, approx (lon, lat)
H_ORBIT = 275.0; RA = RE + H_ORBIT
H_SECO, SECO_DOWNRANGE = 150.0, 15.0                                # km, deg of arc (assumptions, as Flight 14)
DATE, T0_OPEN, T0_CLOSE = '2026-11-13', '22:45:00', '23:59:59'      # 75 min window opens 16:45 CST = 22:45Z, closes 00:00Z

def hms(s):
    h, m, sec = s.split(':'); return int(h)*3600 + int(m)*60 + float(sec)
def fmt(sec):
    s = int(round(sec)); return f'T+{s//3600:02d}:{s%3600//60:02d}:{s%60:02d}'
OFFICIAL = [(hms(l.split('\t')[0]), l.split('\t')[1]) for l in (ROOT/'flight-timeline.txt').read_text().strip().splitlines()]
_KEY = {'Liftoff': 'liftoff', 'Max Q': 'maxq', 'Super Heavy MECO': 'meco', 'Hot-staging': 'hotstage', 'Super Heavy boostback burn start': 'bb0',
        'Super Heavy boostback burn shutdown': 'bb1', 'Super Heavy landing burn start': 'blb', 'Super Heavy landing burn shutdown and catch': 'bcatch',
        'Starship engine cutoff': 'seco', 'Starship orbital insertion burn start': 'ins0', 'Starship orbital insertion burn shutdown': 'ins1',
        'Payload deploy start': 'dep0', 'Payload deploy complete': 'dep1', 'Early Hawaii Return Check': 'chk_haw1', 'Mid Chile Return Check': 'chk_chile',
        'Final Hawaii Return Check': 'chk_haw2', 'Deorbit burn start': 'deorb0', 'Deorbit burn shutdown': 'deorb1', 'Starship entry': 'entry',
        'Starship is transonic': 'transonic', 'Starship is subsonic': 'subsonic', 'Landing burn start': 'lb', 'Landing flip': 'flip',
        'Landing burn 3 to 2 engines': 'lb32', 'Landing burn 2 to 1 engine': 'lb21', 'Starship landing burn shutdown and catch': 'catch'}
T = {}
for t, txt in OFFICIAL:
    for k, v in _KEY.items():
        if txt.startswith(k): T[v] = t
assert len(T) == len(_KEY), set(_KEY.values()) - set(T)
MET_SECO_MIN, MET_INS_MIN = T['seco']/60, T['ins0']/60

# ---- coast ellipse: apogee 275 km; SECO at H_SECO climbing; SECO -> apogee = 18.0 min (Flight 14: 17.3 min)
def _coast_time(rp):
    el = F.Ellipse(RA, rp); nu = el.nu_at_r(RE + H_SECO); return (np.pi - el.M_at_nu(nu))/el.n
RP_COAST = brentq(lambda rp: _coast_time(rp) - (MET_INS_MIN - MET_SECO_MIN), RE - 3000, RE + H_SECO - 1)
COAST = F.Ellipse(RA, RP_COAST); NU_SECO = COAST.nu_at_r(RE + H_SECO); M_SECO = COAST.M_at_nu(NU_SECO)
V_CIRC = np.sqrt(MU/RA)
def geoc2geod(lat_c): return np.degrees(np.arctan(np.tan(np.radians(lat_c))/(1 - EF2)))

class Orbit:
    def __init__(self, inc, dlon, period_scale=0.0):
        """period_scale: the orbital rate is multiplied by (1 + period_scale). flt15_descents solves it with the deorbit delta-v and
        bank so that the modeled catch happens at the official T+22:32:34; it absorbs the uncertain SECO downrange, the mean
        altitude (Flight 14's deorbit telemetry read 270 km against the model's 275) and the Flight 15 ascent differences."""
        self.inc, self.dlon, self.period_scale = inc, dlon, period_scale
        n = np.sqrt(MU/RA**3)*60; k = 1.5*J2*(RE/RA)**2
        self.ci, self.si = np.cos(np.radians(inc)), np.sin(np.radians(inc))
        self.udot = n*(1 + k*(4*self.ci**2 - 1))*(1 + period_scale); self.Odot = -n*k*self.ci; self.T_nodal = 2*np.pi/self.udot   # rad/min, min
        lat_sc = np.arctan((1 - EF2)*np.tan(np.radians(STARBASE[1])))
        self.u0 = np.pi - np.arcsin(np.sin(lat_sc)/self.si)                     # pad, descending branch, at T0
        self.L0 = np.radians(STARBASE[0] + dlon) - np.arctan2(self.ci*np.sin(self.u0), np.cos(self.u0))
        self.u_seco = self.u0 + np.radians(SECO_DOWNRANGE); self.u_ins = self.u_seco + (np.pi - NU_SECO)
        self.P_pad = np.array([np.cos(lat_sc)*np.cos(np.radians(STARBASE[0])), np.cos(lat_sc)*np.sin(np.radians(STARBASE[0])), np.sin(lat_sc)])
        self.P_seco = self._unit(MET_SECO_MIN, self.u_seco)
        d = 1e-3; u_d = self.u_seco + (COAST.nu_at_M(M_SECO + COAST.n*d) - NU_SECO)
        self.V_seco = (self._unit(MET_SECO_MIN + d, u_d) - self.P_seco)/d       # Earth-fixed unit-vector rate at SECO, per minute
        self.v_seco_ground = np.linalg.norm(self.V_seco)/60*RE                  # km/s along the ground
        self.hdot_seco = (COAST.r_at_nu(COAST.nu_at_M(M_SECO + COAST.n*d)) - (RE + H_SECO))/d/60   # km/s
        # ascent ground path: the Flight 14 Hermite curve (pad at rest -> SECO state), tabulated by arc length
        q = np.linspace(0, 1, 6001); h01 = 3*q**2 - 2*q**3; h00 = 1 - h01; h11 = q**3 - q**2
        v = h00[:, None]*self.P_pad + h01[:, None]*self.P_seco + (MET_SECO_MIN*h11)[:, None]*self.V_seco
        v /= np.linalg.norm(v, axis=1)[:, None]; self._asc_pts = v
        self._asc_s = np.r_[0.0, np.cumsum(np.arccos(np.clip((v[1:]*v[:-1]).sum(1), -1, 1)))*RE]
        self.S_total = float(self._asc_s[-1])
        # distance along the path vs time: Hermite through telemetry-like anchors (km, km/s), scaled to end at S_total
        ta = [0, 30, T['maxq'], T['meco'], 200, 300, 400, T['seco']]
        spl = CubicHermiteSpline(ta, [0, 0.3, 3, 55, 165, 500, 1040, 1657], [0, 0.03, 0.13, 1.22, 2.4, 4.4, 6.5, self.v_seco_ground])
        self._scale = self.S_total/float(spl(T['seco'])); self._spl_s = spl
        # height vs time (km, km/s): Max Q 11.5 km, MECO 66 km, then a shallow climb ending on the coast ellipse
        th = [0, 30, T['maxq'], T['meco'], 170, 200, 300, 400, T['seco']]
        self._spl_h = CubicHermiteSpline(th, [0, 2.2, 11.5, 66, 85, 93, 114, 134, H_SECO], [0, 0.16, 0.33, 0.9, 0.35, 0.22, 0.2, 0.19, self.hdot_seco])
        assert (self._spl_s.derivative()(np.linspace(0, T['seco'], 500)) >= -1e-9).all() and (self._spl_h.derivative()(np.linspace(0, T['seco'], 500)) >= -1e-9).all()
    def s_of_t(self, t): return self._spl_s(t)*self._scale                       # km along the ground path
    def v_of_t(self, t): return self._spl_s.derivative()(t)*self._scale         # km/s
    def h_of_t(self, t): return self._spl_h(t)
    def hdot_of_t(self, t): return self._spl_h.derivative()(t)
    def _frame(self, L):
        N = np.array([np.cos(L), np.sin(L), 0.0]); M = np.array([-np.sin(L)*self.ci, np.cos(L)*self.ci, self.si])
        return N, M, np.cross(N, M)
    def _unit(self, met_min, u):
        N, M, _ = self._frame(self.L0 + (self.Odot - OMEGA_E)*met_min); return np.cos(u)*N + np.sin(u)*M
    def lonlat(self, met_min, u):
        L = self.L0 + (self.Odot - OMEGA_E)*np.asarray(met_min, float)
        lon = np.degrees(L + np.arctan2(self.ci*np.sin(u), np.cos(u))); lon = (lon + 180) % 360 - 180
        return lon, geoc2geod(np.degrees(np.arcsin(self.si*np.sin(u))))
    def path_lonlat(self, s_km):
        """lon, lat (geodetic) at ground distance s_km along the ascent path (used by the ascent and the booster)"""
        s = np.clip(np.asarray(s_km, float), 0, self.S_total); idx = np.interp(s, self._asc_s, np.arange(len(self._asc_s)))
        i = np.minimum(idx.astype(int), len(self._asc_s) - 2); w = idx - i
        v = self._asc_pts[i]*(1 - w)[:, None] + self._asc_pts[i + 1]*w[:, None]; v /= np.linalg.norm(v, axis=1)[:, None]
        return (np.degrees(np.arctan2(v[:, 1], v[:, 0])) + 180) % 360 - 180, geoc2geod(np.degrees(np.arcsin(np.clip(v[:, 2], -1, 1))))
    def u_nominal(self, met_min):
        met = np.asarray(met_min, float); u = np.empty_like(met)
        a = met < MET_SECO_MIN; c = (met >= MET_SECO_MIN) & (met < MET_INS_MIN); o = met >= MET_INS_MIN
        u[a] = self.u0 + np.radians(SECO_DOWNRANGE)*self.s_of_t(met[a]*60)/self.S_total
        u[c] = self.u_seco + (COAST.nu_at_M(M_SECO + COAST.n*(met[c] - MET_SECO_MIN)) - NU_SECO)
        u[o] = self.u_ins + self.udot*(met[o] - MET_INS_MIN)
        return u
    def nominal(self, met_min):
        """ground track lon, lat (deg) of the nominal flight up to the deorbit burn"""
        met = np.asarray(met_min, float); lon, lat = self.lonlat(met, self.u_nominal(met)); lon, lat = lon.copy(), lat.copy()
        a = met < MET_SECO_MIN
        if a.any(): lon[a], lat[a] = self.path_lonlat(self.s_of_t(met[a]*60))
        return lon, lat
    def alt_nominal(self, met_min):
        met = np.asarray(met_min, float); h = np.full_like(met, H_ORBIT)
        a = met < MET_SECO_MIN; c = (met >= MET_SECO_MIN) & (met < MET_INS_MIN)
        h[a] = self.h_of_t(met[a]*60); h[c] = COAST.r_at_nu(COAST.nu_at_M(M_SECO + COAST.n*(met[c] - MET_SECO_MIN))) - RE
        return h
    def orbit_state(self, met_s):
        """ECEF position (km) and Earth-relative velocity (km/s) on the circular orbit at MET met_s (after insertion)"""
        m = met_s/60; u = float(self.u_nominal(np.array([m]))[0]); N, M, _ = self._frame(self.L0 + (self.Odot - OMEGA_E)*m)
        r = RA*(np.cos(u)*N + np.sin(u)*M); v_in = V_CIRC*(-np.sin(u)*N + np.cos(u)*M)
        return r, v_in - np.cross([0, 0, OMEGA_E/60], r)
    def seco_state(self):
        """ECEF state on the coast ellipse at SECO (for the no-insertion contingency)"""
        u = self.u_seco; N, M, _ = self._frame(self.L0 + (self.Odot - OMEGA_E)*MET_SECO_MIN)
        rh, th = np.cos(u)*N + np.sin(u)*M, -np.sin(u)*N + np.cos(u)*M
        r0 = (RE + H_SECO)*rh; c = np.sqrt(MU/COAST.p); nu = NU_SECO
        v_in = c*COAST.e*np.sin(nu)*rh + c*(1 + COAST.e*np.cos(nu))*th
        return T['seco'], r0, v_in - np.cross([0, 0, OMEGA_E/60], r0)

import json as _json
PHASE_FILE = ROOT/'work'/'phase.json'
PERIOD_SCALE = _json.load(open(PHASE_FILE))['period_scale'] if PHASE_FILE.exists() else 0.0
orbit = Orbit(F.inc_fit, F.dlon_fit, PERIOD_SCALE)

# ---------------- Super Heavy: MECO -> boostback -> fall -> landing burn -> catch ----------------
BETA_BOOSTER = 3400.0      # kg/m^2: ~270 t falling engines-first, C_D A ~ 77 m^2 with the grid fins (assumption)
H_LB = 4.5                 # km, landing-burn start height; the start point east of the tower follows from the 25 s burn to rest
S_LB = 0.0
def _booster_rhs(t, y, a_bb, th_bb):
    s, h, vs, vh = y; g = MU/(RE + h)**2; rho = atmo(h)[0]; V = np.hypot(vs, vh)
    D = 0.5*rho*(V*1000)**2/BETA_BOOSTER/1000 if V > 0 else 0.0                # km/s^2
    acs, ach = (-D*vs/V, -D*vh/V) if V > 0 else (0.0, 0.0)
    if T['bb0'] <= t < T['bb1']: acs -= a_bb*np.cos(th_bb); ach += a_bb*np.sin(th_bb)
    return [vs, vh, acs, ach - g]
def _booster_fall(p, record=False):
    a_bb, th_bb = p[0]/1000, np.radians(p[1])
    y0 = [float(orbit.s_of_t(T['meco'])), float(orbit.h_of_t(T['meco'])), float(orbit.v_of_t(T['meco'])), float(orbit.hdot_of_t(T['meco']))]
    sol = solve_ivp(_booster_rhs, (T['meco'], T['blb']), y0, args=(a_bb, th_bb), method='DOP853', rtol=1e-9, atol=1e-9, max_step=2.0, dense_output=record)
    return sol
def solve_booster():
    """boostback acceleration (m/s^2) and pitch (deg) so the fall reaches H_LB at S_LB at T+06:10; S_LB is iterated so the
    25 s landing burn, a linear deceleration of the arrival velocity, ends at rest on the tower (the booster arrives from the Gulf)"""
    global S_LB
    for _ in range(3):
        best = None
        for x0 in ([40.0, 10.0], [30.0, 30.0]):
            r = least_squares(lambda p: (lambda y: [y[0] - S_LB, y[1] - H_LB])(_booster_fall(p).y[:, -1]), x0, bounds=([5, -30], [120, 80]), x_scale=[10, 10], diff_step=1e-3)
            if best is None or r.cost < best.cost: best = r
        vs_lb = _booster_fall(best.x).y[2, -1]; S_LB = -vs_lb*(T['bcatch'] - T['blb'])/2
    return best.x
def booster_track(p, dt=1.0):
    """1 s track (t, s, h, vs, vh) from MECO to the catch; the landing burn is a cubic Hermite to rest at the tower"""
    sol = _booster_fall(p, record=True)
    t1 = np.arange(T['meco'], T['blb'] + 1e-9, dt); y1 = sol.sol(t1)
    s_lb, h_lb, vs_lb, vh_lb = sol.y[:, -1]
    lb = CubicHermiteSpline([T['blb'], T['bcatch']], [[s_lb, h_lb], [0.0, 0.0]], [[vs_lb, vh_lb], [0.0, 0.0]])    # linear deceleration to rest
    t2 = np.arange(T['blb'] + dt, T['bcatch'] + 1e-9, dt); y2 = lb(t2).T; v2 = lb.derivative()(t2).T
    t = np.r_[t1, t2]; s = np.r_[y1[0], y2[0]]; h = np.r_[y1[1], np.maximum(y2[1], 0)]; vs = np.r_[y1[2], v2[0]]; vh = np.r_[y1[3], v2[1]]
    i_ap = int(np.argmax(h)); i_max = int(np.argmax(s))
    return dict(t=t, s=s, h=h, vs=vs, vh=vh, apogee=(t[i_ap], h[i_ap]), max_downrange=(t[i_max], s[i_max]), v_lb=np.hypot(vs_lb, vh_lb),
                dv_bb=p[0]*(T['bb1'] - T['bb0']), a_bb=p[0], th_bb=p[1], s_lb=s_lb, sol=sol)

if __name__ == '__main__':
    o = orbit
    print(f'plane i={o.inc:.3f} node offset {o.dlon:+.3f} deg; nodal period {o.T_nodal:.3f} min (period scale {o.period_scale:+.5f}); shift {np.degrees((OMEGA_E - o.Odot)*o.T_nodal):.3f} deg W/rev')
    print(f'coast ellipse: perigee {RP_COAST - RE:.0f} km, e={COAST.e:.4f}; SECO ground speed {o.v_seco_ground:.3f} km/s, climb {o.hdot_seco*1000:.0f} m/s; path length {o.S_total:.0f} km')
    for t in [0, 30, T['maxq'], 100, T['meco'], T['hotstage'], 200, 300, 400, T['seco']]:
        lo, la = o.path_lonlat([o.s_of_t(t)]); print(f'  {fmt(t)}  s {o.s_of_t(t):7.1f} km  h {o.h_of_t(t):6.1f} km  v_ground {o.v_of_t(t):5.2f} km/s  climb {o.hdot_of_t(t)*1000:5.0f} m/s  {la[0]:.3f},{lo[0]:.3f}')
    for k in ['seco', 'ins0', 'dep0', 'dep1', 'chk_haw1', 'chk_chile', 'chk_haw2', 'deorb0', 'entry', 'catch']:
        lo, la = o.nominal(np.array([T[k]/60])); u = np.degrees(o.u_nominal(np.array([T[k]/60]))[0]) % 360
        print(f'{k:9s} {fmt(T[k])}  u {u:7.2f}  {la[0]:7.2f}, {lo[0]:8.2f}  h {o.alt_nominal(np.array([T[k]/60]))[0]:.1f}')
    asc = []
    mm = np.arange(0, 24*60, 1/6); lo, la = o.nominal(mm); i = np.where(np.diff(np.sign(la)) > 0)[0]
    print('ascending nodes:', [(fmt(mm[j]*60), round(float(lo[j]), 2)) for j in i])
    p = solve_booster(); b = booster_track(p)
    print(f'booster: boostback a={b["a_bb"]:.1f} m/s^2 pitch {b["th_bb"]:.1f} deg, dv {b["dv_bb"]:.0f} m/s; apogee {b["apogee"][1]:.1f} km at {fmt(b["apogee"][0])}; '
          f'max downrange {b["max_downrange"][1]:.1f} km at {fmt(b["max_downrange"][0])}; landing-burn start {b["v_lb"]*1000:.0f} m/s at {H_LB} km, {b["s_lb"]:.1f} km east of the tower')
    for t in [T['meco'], T['bb0'], T['bb1'], 200, 233, 300, T['blb'], 385, T['bcatch']]:
        j = int(np.argmin(np.abs(b['t'] - t))); print(f'  {fmt(t)}  s {b["s"][j]:6.1f} km  h {b["h"][j]:6.1f} km  v ({b["vs"][j]*1000:6.0f}, {b["vh"][j]*1000:6.0f}) m/s')
