"""Starship Flight 15 descents on the Flight 14-calibrated ship (ift14_entry: 3-DOF point mass, rotating WGS-84 Earth, J2, US76
atmosphere, Newtonian belly-first aerodynamics with the Flight 14 angle of attack and ballistic coefficient, bank steering, a
10 s one-Raptor deorbit burn), every one started from the Flight 15 orbit (flt15_orbit):

  nominal   deorbit burn T+21:35:57 (official); landing burn 18 s (official). Two stages: (1) delta-v, a constant bank and the orbit
            period scale solved for the Pad 2 tower at the official catch time (this fixes the orbit for every descent); (2) with that
            orbit, delta-v and a two-stage bank solved so that the path runs along the entry corridor over Mexico (NOTAM B4288/26,
            areas C2 and C3) and covers the routed distance to the tower by the official catch time
  hawaii1   early Hawaii: orbit-2 burn at T+02:12:42 (Flight 14's T+02:12:00 re-timed by the 42 s later insertion), solved for
            Flight 14's splashdown point north of Hawaii (North Pacific zone)
  chile     mid-mission: orbit-7 burn at T+08:53:00 (Flight 14's planned T+08:52:18 re-timed), delta-v and bank solved for
            Flight 14's planned landing point in the South Pacific zone west of Chile
  hawaii2   late Hawaii: orbit 16, the pass after the Starbase pass (HYDROPAC 3256/26 opens 23:32 after the earliest launch, one
            orbit after the nominal entry): a Flight 14-like 69 m/s burn about one orbit after the nominal one, burn time and bank
            solved for Flight 14's splashdown point north-east of Hawaii, inside the area's eastern lobe
  indian    no insertion burn: ballistic coast from SECO with the payload aboard, bank steering to the middle of the Indian Ocean zone

The nominal track east of 103.5 W is re-shaped kinematically, at equal length, to the published areas and SpaceX's description of
the approach: through the east end of the Mexican corridor on the Rio Grande, over the Laguna Madre 57 km north of the South Padre
Island town (ship-catch area NAVAREA IV 1081/26: Laguna Madre and Padre Island from 26.8 N down to the pad, the Gulf out to 96.4 W),
out over the Gulf, final turn south-west to the tower; heights, speeds and times stay those of the simulation.
Results are cached in descents.json; run with --solve to recompute (about 15 minutes), or with --corridor to redo only the nominal's
corridor stage (about 4 minutes; the orbit and the contingencies stay)."""
import sys, json, numpy as np
from scipy.optimize import least_squares, brentq
from scipy.interpolate import CubicHermiteSpline
from shapely.geometry import LineString, Point
from shapely.ops import unary_union
sys.path.insert(0, '/Users/mickey/sda/starship-ITF-14/work'); sys.path.insert(0, '/Users/mickey/launch/Starship-FLT15/work')
import ift14_fit as F, ift14_entry as E
from ift14_atmo import atmo
import flt15_orbit as O, flt15_zones as Z
T, fmt, orbit = O.T, O.fmt, O.orbit
CACHE = O.ROOT/'work'/'descents.json'
X14 = json.load(open(F.ROOT/'work'/'entry_solution.json'))['x']; DV14, SIG14, _, BETA, ALPHA = X14
S_REF = E.S_REF; M_ENTRY = BETA*(4/3*np.sin(np.radians(ALPHA))**3)*S_REF/1000; PAYLOAD_T = 50.0
BURN_S, LB_CATCH_S, LB_SEA_S = 10.0, T['catch'] - T['lb'], 19.0
TOWER = (O.TOWER[1], O.TOWER[0])                       # lat, lon
HAWAII14 = (25.499295, -155.427536)                     # Flight 14 splashdown (reported)
DV_LATE = 0.069                                         # km/s, late-Hawaii deorbit burn: the Flight 14 / early-Hawaii value
T_BURN_HAW1 = 2*3600 + 12*60 + 42.0
T_BURN_CHILE = 8*3600 + 53*60 + 0.0
PLACES = {'Lihue (Kauai)': (21.98, -159.37), 'Honolulu': (21.31, -157.86), 'Hilo': (19.72, -155.09), 'Hermosillo': (29.07, -110.96),
          'Chihuahua': (28.63, -106.08), 'Ciudad Juarez / El Paso': (31.74, -106.49), 'Monterrey': (25.67, -100.31), 'Laredo': (27.51, -99.51),
          'Monclova': (26.9, -101.42), 'Reynosa / McAllen': (26.1, -98.3), 'Harlingen': (26.19, -97.70), 'Brownsville': (25.9, -97.5),
          'South Padre Island': (26.11, -97.17), 'Corpus Christi': (27.8, -97.4), 'Raymondville': (26.48, -97.78), 'Port Mansfield': (26.55, -97.43),
          'Piedras Negras / Eagle Pass': (28.70, -100.52), 'Nuevo Laredo': (27.49, -99.52), 'Ojinaga / Presidio': (29.56, -104.41)}

class Sim15(E.Sim):
    """Flight 14's ship started from a given ECEF state; optional second bank sigma2 after t_rev (default: a reversal)"""
    def __init__(self, dv, sigma, t_rev, t_burn, start, sigma2=None, lb_dur=LB_SEA_S, beta=BETA):
        super().__init__(dv, sigma, t_rev, beta, ALPHA, t_burn=t_burn, start=start)
        self.sigma2 = None if sigma2 is None else np.radians(sigma2); self.lb_dur = lb_dur
    def bank(self, t, M):
        s = self.sigma if t < self.t_rev else (-self.sigma if self.sigma2 is None else self.sigma2)
        return s*np.clip((M - 1.5)/1.5, 0, 1)
    def run(self, **kw):
        E.BURN_DUR = BURN_S; E.LB_DUR = self.lb_dur; return super().run(**kw)

def orbit_start(t_burn, eps=None):
    ob = orbit if eps is None else O.Orbit(F.inc_fit, F.dlon_fit, eps); return (t_burn, *ob.orbit_state(t_burn))
def miss_km(s, target):
    la, lo = s.splash; return (la - target[0])*111.2, ((lo - target[1] + 540) % 360 - 180)*111.2*np.cos(np.radians(la))
def nan_lb(s): return np.isnan(getattr(s, 't_lb', np.nan))

def fit_point(make, x0s, target, bounds, x_scale):
    best = None
    for x0 in x0s:
        def res(p):
            s = make(p).run()
            if nan_lb(s): return [1e3, 1e3]
            n, e = miss_km(s, target); return [n/5, e/5]
        r = least_squares(res, x0, x_scale=x_scale, bounds=bounds, diff_step=2e-3, max_nfev=120)
        print(f'    start {np.round(x0, 3)} -> {np.round(r.x, 4)} miss {np.round(r.fun*5, 2)} km', flush=True)
        if best is None or r.cost < best.cost: best = r
    return best.x

def solve_nominal():
    """delta-v, bank and the orbit period scale for the tower at the official catch time"""
    tb = T['deorb0']; best = None
    def res(p):
        s = Sim15(p[0], p[1], 1e9, tb, orbit_start(tb, p[2]), lb_dur=LB_CATCH_S).run()
        if nan_lb(s): return [1e3, 1e3, 1e3]
        n, e = miss_km(s, TOWER); return [n/5, e/5, (s.t_splash - T['catch'])/5]
    for x0 in ([0.0666, 19.5, 0.0014], [0.06, 25.0, 0.001]):
        r = least_squares(res, x0, x_scale=[0.005, 10, 5e-4], bounds=([0.02, -89, -0.01], [0.15, 89, 0.01]), diff_step=2e-3, max_nfev=150)
        print(f'    start {np.round(x0, 4)} -> {np.round(r.x, 5)} resid [N km, E km, s] {np.round(r.fun*5, 2)}', flush=True)
        if best is None or r.cost < best.cost: best = r
    x = best.x
    return dict(dv=float(x[0]), sigma=float(x[1]), t_rev=1e9, t_burn=tb, lb_dur=LB_CATCH_S, period_scale=float(x[2]))
def solve_hawaii1():
    st = orbit_start(T_BURN_HAW1)
    x = fit_point(lambda p: Sim15(p[0], p[1], 1e9, T_BURN_HAW1, st), ([DV14, SIG14], [0.06, 25.0]), HAWAII14, ([0.02, -89], [0.15, 89]), [0.005, 10])
    return dict(dv=float(x[0]), sigma=float(x[1]), t_rev=1e9, t_burn=T_BURN_HAW1, lb_dur=LB_SEA_S)
CHILE14 = (-29.9, -86.9)                                 # Flight 14's planned landing point, inside the zone's eastern part
def solve_chile():
    st = orbit_start(T_BURN_CHILE)
    x = fit_point(lambda p: Sim15(p[0], p[1], 1e9, T_BURN_CHILE, st), ([0.069, 0.0], [0.065, -15.0]), CHILE14, ([0.02, -89], [0.15, 89]), [0.005, 10])
    return dict(dv=float(x[0]), sigma=float(x[1]), t_rev=1e9, t_burn=T_BURN_CHILE, lb_dur=LB_SEA_S)
def solve_hawaii2():
    """late Hawaii on orbit 16, the pass after the Starbase pass: HYDROPAC 3256/26 runs from Micronesia to north-east of Hawaii and
    its window opens one orbit after the nominal entry, so the burn is the Flight 14-like 69 m/s burn about one nodal period after
    the nominal one; burn time and bank solved for Flight 14's splashdown point, which lies in the area's eastern lobe (the Final
    Hawaii Return Check at T+20:41:55 and the nominal deorbit opportunity precede it)"""
    best = None; tb0 = T['deorb0'] + orbit.T_nodal*60
    def res(p):
        tb = float(p[0]); s = Sim15(DV_LATE, p[1], 1e9, tb, orbit_start(tb)).run()
        if nan_lb(s): return [1e3, 1e3]
        n, e = miss_km(s, HAWAII14); return [n/5, e/5]
    for x0 in ([tb0 - 6*60, 20.0], [tb0 - 4*60, 25.0]):
        r = least_squares(res, x0, x_scale=[60, 10], bounds=([tb0 - 20*60, -89], [tb0 + 5*60, 89]), diff_step=2e-3, max_nfev=120)
        print(f'    start {np.round(x0, 3)} -> {np.round(r.x, 4)} miss {np.round(r.fun*5, 2)} km', flush=True)
        if best is None or r.cost < best.cost: best = r
    return dict(dv=DV_LATE, sigma=float(best.x[1]), t_rev=1e9, t_burn=float(best.x[0]), lb_dur=LB_SEA_S)

def solve_indian():
    cl = Z.centerline(Z.ZD['indian']); lon_t = 0.5*(cl[0, 0] + cl[-1, 0]); lat_t = float(np.interp(lon_t, cl[:, 0], cl[:, 1]))
    st = orbit.seco_state(); beta = BETA*(M_ENTRY + PAYLOAD_T)/M_ENTRY
    x = fit_point(lambda p: Sim15(0.0, p[0], p[1], 1e9, st, beta=beta), ([60.0, 2900.0], [75.0, 3100.0]), (lat_t, lon_t), ([0, 2000], [89, 4500]), [10, 200])
    return dict(dv=0.0, sigma=float(x[0]), t_rev=float(x[1]), t_burn=1e9, lb_dur=LB_SEA_S, beta=float(beta), target=[lat_t, lon_t])

def make_sim(p, start):
    return Sim15(p['dv'], p['sigma'], p['t_rev'], p['t_burn'], start, sigma2=p.get('sigma2'), lb_dur=p['lb_dur'], beta=p.get('beta', BETA))

def sample(sim, t_from, dt=2.0):
    """dt track from t_from to the landing (last sample = the landing point), events, loads and the lon/lat/alt/speed arrays"""
    s = sim.run(record=True); sol = s.sol; LB = sim.lb_dur
    def state(t):
        if t <= s.t_lb: y = sol.sol(t); return y[:3], y[3:]
        y = sol.sol(s.t_lb); r, v = y[:3], y[3:]; up = r/np.linalg.norm(r)
        vz = -(v @ up); vh = v - (v @ up)*up; tau = min(t - s.t_lb, LB)
        h = max(E.geodetic(r)[2] - vz*tau + vz*tau*tau/(2*LB), 0.0)
        la, lo, _ = E.geodetic(r + vh*(tau - tau*tau/(2*LB))); return E.ecef(la, lo, h), (vh - vz*up)*max(1 - tau/LB, 0)
    grid = np.arange(t_from, s.t_splash - 0.5*dt, dt)
    lat, lon, alt, spd = [], [], [], []
    for t in grid:
        r, v = state(t); la, lo, h = E.geodetic(r); lat.append(la); lon.append(lo); alt.append(max(h, 0.0)); spd.append(np.linalg.norm(v))
    lat.append(s.splash[0]); lon.append(s.splash[1]); alt.append(0.0); spd.append(0.0)
    met = np.r_[grid, s.t_splash]
    tt = np.arange(max(t_from, s.t_ei - 600), s.t_lb, 1.0); M, Q, G, H = [], [], [], []
    for t in tt:
        y = sol.sol(t); h = E.geodetic(y[:3])[2]; rho, Tk, c = atmo(h); V = np.linalg.norm(y[3:])*1000
        ld, cdr = E.aero(V/c, ALPHA); D = 0.5*rho*V*V*cdr/sim.beta if h < 200 else 0.0
        M.append(V/c); Q.append(np.sqrt(rho)*V**3); G.append(D*np.sqrt(1 + ld*ld)/9.80665); H.append(h)
    M, Q, G, H = map(np.array, (M, Q, G, H))
    below = lambda lev: float(tt[np.argmax((M < lev) & (tt > s.t_ei))])
    alt_at = lambda t: float(E.geodetic(sol.sol(min(max(t, tt[0]), s.t_lb))[:3])[2])
    ei = E.geodetic(sol.sol(s.t_ei)[:3]); bp = E.geodetic(sol.sol(min(sim.t_burn, s.t_lb))[:3]) if sim.t_burn < 1e8 else E.geodetic(sim.start[1])
    return dict(met=met.tolist(), lon=lon, lat=lat, alt=alt, speed=spd, dt=dt, t_burn=float(sim.t_burn) if sim.t_burn < 1e8 else None,
                t_ei=float(s.t_ei), ei=[ei[1], ei[0]], burn=[bp[1], bp[0]], splash=[float(s.splash[1]), float(s.splash[0])], t_splash=float(s.t_splash),
                peak_heating=float(tt[int(np.argmax(Q))]), peak_heating_alt=round(float(H[int(np.argmax(Q))]), 1), max_g=round(float(G.max()), 2),
                t_max_g=float(tt[int(np.argmax(G))]), transonic=below(1.2), subsonic=below(1.0), mach5=below(5.0), mach3=below(3.0),
                landing_burn=float(s.t_lb), alt_at_official_entry=alt_at(T['entry']) if sim.t_burn == T['deorb0'] else None, near=nearby(met, lon, lat, alt))
def nearby(met, lon, lat, alt):
    """closest approach (km, MET, height) of a ground track to the places of interest within 400 km"""
    near = {}
    for name, (pla, plo) in PLACES.items():
        d = F.gc_dist_deg(np.array(lon), np.array(lat), plo, pla)*111.2; i = int(np.argmin(d))
        if d[i] < 400: near[name] = dict(km=round(float(d[i]), 1), met=float(met[i]), alt=round(float(alt[i]), 1))
    return near

# ---- the entry corridor over Mexico and the final approach: the published areas and SpaceX's description laid over the simulated profile ----
CORRIDOR = unary_union([Z.ZD['mexC2']['poly'], Z.ZD['mexC3']['poly']])     # NOTAM B4288/26: Pacific off Baja California (120.6 W) to the Rio Grande (99.9 W)
CORR_LONS = np.arange(-119.5, -100.4, 1.0)
def _cut(L): c = CORRIDOR.intersection(LineString([(L, 0), (L, 60)])); return c.bounds[1], c.bounds[3]
CORR_EDGE = np.array([_cut(L) for L in CORR_LONS]); CORR_CEN = CORR_EDGE.mean(axis=1); CORR_HALF = (CORR_EDGE[:, 1] - CORR_EDGE[:, 0])*55.6   # km
W_BORDER, AZ_BORDER = (27.942, -99.942), 108.0        # middle of the corridor's east end on the Rio Grande (28-04N 100-01W to 27-49N 099-52W), 64 km north-west of Laredo; its axis there
W_LAGUNA, W_GULF = (26.62, -97.35), (26.18, -96.85)   # over the Laguna Madre and Padre Island 57 km north of the South Padre Island town (the hazard
START_LON = -103.5                                     # area's northern tip is 26.77 N 97.47 W); over the Gulf 36 km north-east of the tower before the final turn
def path_km(lon, lat):
    lon, lat = np.asarray(lon), np.asarray(lat); dl = (np.diff(lon) + 540) % 360 - 180
    return np.r_[0, np.cumsum(np.hypot(dl*111.2*np.cos(np.radians(0.5*(lat[1:] + lat[:-1]))), np.diff(lat)*111.2))]
def route(lon0, lat0, lon1, lat1):
    """the re-shaped ground path from a point of the simulated track (heading toward lon1, lat1): corridor east end, Laguna Madre,
    Gulf, tower, on a chord-length Hermite with unit tangents -> dense lon, lat and the cumulative length in km"""
    lat_t, lon_t = TOWER; kx, ky = 111.2*np.cos(np.radians(lat_t)), 111.2; xy = lambda la, lo: np.array([(lo - lon_t)*kx, (la - lat_t)*ky])
    p0 = xy(lat0, lon0); d0 = xy(lat1, lon1) - p0; d0 /= np.linalg.norm(d0)
    P = np.array([p0, xy(*W_BORDER), xy(*W_LAGUNA), xy(*W_GULF), np.zeros(2)]); u = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    az = lambda a: np.array([np.sin(np.radians(a)), np.cos(np.radians(a))])
    Cc = CubicHermiteSpline(u, P, np.array([d0, az(AZ_BORDER), az(140), az(200), az(235)]))(np.linspace(0, u[-1], 8000))   # u is the chord length, so |dP/du| ~ 1
    lo, la = lon_t + Cc[:, 0]/kx, lat_t + Cc[:, 1]/ky; return lo, la, path_km(lo, la)
def reshape(met, lon, lat, t_from):
    """lay the route over a sampled track from its first point east of START_LON (after t_from): new lon, lat, the start index, and
    the route length over the simulated distance left at that point (1 when the simulation flies the routed distance)"""
    met, lon, lat = np.asarray(met), np.array(lon, float), np.array(lat, float)
    k = int(np.argmax((met > t_from) & (lon > START_LON) & (lon < -60))); dist = path_km(lon, lat); rem = dist[-1] - dist[k]
    rl, ra, L = route(lon[k], lat[k], lon[k + 1], lat[k + 1]); stretch = L[-1]/rem
    d = np.minimum((dist[k:] - dist[k])*stretch, L[-1]); lon[k:] = np.interp(d, L, rl); lat[k:] = np.interp(d, L, ra)
    return lon, lat, k, float(stretch), float(rem)
def approach_overlay(tr):
    """re-shape the ground path of the nominal descent east of START_LON so that it leaves Mexico through the middle of the corridor's
    east end, reaches the Laguna Madre well north of South Padre Island heading south-east, runs out over the Gulf and turns south-west
    to the tower, inside the ship-catch hazard area; heights, speeds and times stay those of the simulation. The simulation is solved to
    fly this routed distance (solve_corridor), so the horizontal spacing is kept (overlay_stretch = 1 to the solver's tolerance)"""
    lon, lat, k, stretch, rem = reshape(tr['met'], tr['lon'], tr['lat'], tr['t_ei'])
    tr['sim_end'] = [float(tr['lon'][-1]), float(tr['lat'][-1])]                # where the un-routed simulation comes down
    tr['lon'], tr['lat'] = lon.tolist(), lat.tolist(); tr['overlay_from'] = float(tr['met'][k]); tr['overlay_km'] = rem; tr['overlay_stretch'] = stretch; tr['splash'] = [float(lon[-1]), float(lat[-1])]
    tr['near'] = nearby(np.array(tr['met']), lon, lat, np.array(tr['alt']))     # closest approaches of the re-shaped path
    return tr
def corridor_offsets(lon, lat):
    """north-south offset (km, + north) of a track from the corridor's centreline at CORR_LONS"""
    lon, lat = np.asarray(lon), np.asarray(lat); m = (lon > -125) & (lon < -98); o = np.argsort(lon[m]); return (np.interp(CORR_LONS, lon[m][o], lat[m][o]) - CORR_CEN)*111.2
def solve_corridor(direct):
    """second stage of the nominal: the orbit (period scale) stays as solved for the direct path; delta-v, a first bank, a second bank
    and the switch time are solved so that (1) the catch is at the official time, (2) the simulated distance east of START_LON equals the
    routed distance to the tower (corridor east end, Laguna Madre, Gulf), and (3) the path follows the centreline of the Mexican entry
    corridor, each offset weighted by the local half-width. A single bank cannot do this: the corridor keeps within ~30 km of the orbital
    ground track until 120 W and has turned 230 km south of it by the Rio Grande"""
    tb = T['deorb0']; st = orbit_start(tb); best = None; direct = direct.get('direct', direct)
    def run(p): return Sim15(p[0], p[1], p[3], tb, st, sigma2=p[2], lb_dur=LB_CATCH_S).run(record=True)
    def track(s):
        tt = np.arange(s.t_ei, s.t_lb, 2.0); ys = s.sol.sol(tt); g = np.array([E.geodetic(ys[:3, i]) for i in range(len(tt))])
        return np.r_[tt, s.t_splash], np.r_[g[:, 1], s.splash[1]], np.r_[g[:, 0], s.splash[0]]
    def res(p):
        s = run(p)
        if nan_lb(s): return np.r_[1e3, 1e3, np.full(len(CORR_LONS), 1e2)]
        tt, lo, la = track(s); _, _, _, stretch, rem = reshape(tt, lo, la, s.t_ei)
        return np.r_[(s.t_splash - T['catch'])/0.05, (stretch - 1)*rem/0.2, corridor_offsets(lo[:-1], la[:-1])/(0.3*CORR_HALF)]
    for x0 in ([direct['dv'] + 1.7e-4, 11.0, 27.0, T['entry'] + 500], [direct['dv'], 7.0, 26.0, T['entry'] + 390]):
        r = least_squares(res, x0, x_scale=[0.003, 8, 8, 120], bounds=([0.03, -60, -60, T['entry']], [0.12, 60, 80, T['catch'] - 120]), diff_step=[2e-3, 5e-3, 5e-3, 3e-4], max_nfev=80)
        print(f'    start {np.round(x0, 3)} -> dv {r.x[0]*1000:.3f} m/s bank {r.x[1]:.2f} -> {r.x[2]:.2f} at {fmt(r.x[3])}; catch {r.fun[0]*0.05:+.3f} s, routed distance {r.fun[1]*0.2:+.3f} km, '
              f'corridor offsets {np.round(r.fun[2:]*0.3*CORR_HALF).astype(int)} km, cost {r.cost:.2f}', flush=True)
        if best is None or r.cost < best.cost: best = r
    x = best.x
    return dict(dv=float(x[0]), sigma=float(x[1]), sigma2=float(x[2]), t_rev=float(x[3]), direct=dict(dv=direct['dv'], sigma=direct['sigma']))

def build(solve=False, corridor=False):
    cache = json.load(open(CACHE)) if CACHE.exists() and not solve else {}
    sol = cache.get('solutions', {})
    if corridor and 'nominal' in sol: sol['nominal'].pop('sigma2', None)      # redo only the nominal's corridor stage
    global orbit
    for key, fn in [('nominal', solve_nominal), ('hawaii1', solve_hawaii1), ('chile', solve_chile), ('hawaii2', solve_hawaii2), ('indian', solve_indian)]:
        if key not in sol:
            print(f'solving {key} ...', flush=True); sol[key] = fn(); print(f'  {key}: {sol[key]}', flush=True)
            cache['solutions'] = sol; json.dump(cache, open(CACHE, 'w'), default=float)
        if key == 'nominal':                     # every other descent starts from the orbit with the solved period scale
            json.dump(dict(period_scale=sol['nominal']['period_scale']), open(O.PHASE_FILE, 'w'))
            O.PERIOD_SCALE = sol['nominal']['period_scale']; O.orbit = orbit = O.Orbit(F.inc_fit, F.dlon_fit, O.PERIOD_SCALE)
            if 'sigma2' not in sol['nominal']:   # second stage, on the orbit just fixed: along the Mexican corridor and the routed approach
                print('solving nominal, corridor stage ...', flush=True); sol['nominal'].update(solve_corridor(sol['nominal'])); print(f"  nominal: {sol['nominal']}", flush=True)
                cache.get('tracks', {}).pop('nominal', None); cache['solutions'] = sol; json.dump(cache, open(CACHE, 'w'), default=float)
    tr = {} if solve else cache.get('tracks', {})
    if any(key not in tr for key in sol):
        for key, p in sol.items():
            if key in tr: continue
            start = orbit.seco_state() if key == 'indian' else orbit_start(p['t_burn'])
            tr[key] = sample(make_sim(p, start), T['seco'] if key == 'indian' else p['t_burn'])
            print(f'  sampled {key}: {len(tr[key]["met"])} pts, splash {fmt(tr[key]["t_splash"])}', flush=True)
            if key == 'nominal': tr[key] = approach_overlay(tr[key])
        cache['tracks'] = tr; cache['vehicle'] = dict(alpha=ALPHA, beta=BETA, m_entry_t=M_ENTRY, payload_t=PAYLOAD_T, burn_s=BURN_S)
        json.dump(cache, open(CACHE, 'w'), default=float)
    return cache

if __name__ == '__main__':
    C = build(solve='--solve' in sys.argv, corridor='--corridor' in sys.argv)
    ZONE = {'nominal': 'catch', 'hawaii1': 'npac1', 'hawaii2': 'late', 'chile': 'chile', 'indian': 'indian'}
    for key, tr in C['tracks'].items():
        p = C['solutions'][key]; zk = ZONE[key]
        inz = Z.contains(Z.ZD[zk], *tr['splash'])
        print(f"{key:8s} dv {p['dv']*1000:5.1f} m/s bank {p['sigma']:6.1f} deg{' then %.1f deg from %s' % (p['sigma2'], fmt(p['t_rev'])) if p.get('sigma2') is not None else ''}{' period scale %+.5f (T_nodal %.3f min)' % (p['period_scale'], orbit.T_nodal) if 'period_scale' in p else ''} rev {fmt(p['t_rev']) if p['t_rev'] < 1e8 and p.get('sigma2') is None else '-'} | burn {fmt(tr['t_burn']) if tr['t_burn'] else '-'} @ {tr['burn'][1]:.1f},{tr['burn'][0]:.1f} "
              f"| 120 km {fmt(tr['t_ei'])} @ {tr['ei'][1]:.1f},{tr['ei'][0]:.1f} | peak heat {fmt(tr['peak_heating'])} {tr['peak_heating_alt']} km | max {tr['max_g']} g {fmt(tr['t_max_g'])} "
              f"| M5 {fmt(tr['mach5'])} transonic {fmt(tr['transonic'])} subsonic {fmt(tr['subsonic'])} LB {fmt(tr['landing_burn'])} | end {fmt(tr['t_splash'])} @ {tr['splash'][1]:.3f},{tr['splash'][0]:.3f} in zone {inz}")
        if tr.get('alt_at_official_entry') is not None: print(f"          altitude at the official entry {fmt(T['entry'])}: {tr['alt_at_official_entry']:.1f} km")
        cr = Z.crossings(Z.ZD[zk], tr['met'], tr['lon'], tr['lat']); (a, b), (c, d) = Z.window_met(Z.ZD[zk])
        if cr:
            j = int(np.searchsorted(tr['met'], cr[0][0])); print(f"          in {Z.ZD[zk]['name']} from {fmt(cr[0][0])} ({tr['alt'][j]:.0f} km) to {fmt(cr[-1][1])}; "
                                                               f"window {Z.ZD[zk]['t0']}-{Z.ZD[zk]['t1']}Z = {fmt(a)}-{fmt(b)} after a 22:45Z launch, {fmt(c)}-{fmt(d)} after 00:00Z")
        else: print(f"          never inside {Z.ZD[zk]['name']}")
        for n, q in tr['near'].items(): print(f"          near {n:24s} {q['km']:6.1f} km at {fmt(q['met'])}  alt {q['alt']:.1f} km")
    n = C['tracks']['nominal']; met, lon, lat, alt = map(np.array, (n['met'], n['lon'], n['lat'], n['alt']))
    for zk in ('mexC2', 'mexC3'):
        z = Z.ZD[zk]; (a, b), (c, d) = Z.window_met(z)
        for t_in, t_out in Z.crossings(z, met, lon, lat):
            i, j = np.searchsorted(met, [t_in, t_out]); print(f"nominal in {z['name']} from {fmt(t_in)} ({lon[i]:.2f}, {lat[i]:.2f}, {alt[i]:.0f} km) to {fmt(t_out)} ({lon[j]:.2f}, {lat[j]:.2f}, {alt[j]:.0f} km); "
                                                              f"window {z['t0']}-{z['t1']}Z = {fmt(a)}-{fmt(b)} after a 22:45Z launch, {fmt(c)}-{fmt(d)} after 00:00Z")
    m = (met > n['t_ei']) & (lon > -121) & (lon < -99.9); off = corridor_offsets(lon[m], lat[m])
    print('  offset from the corridor centreline (km, + north) / half-width, every 1 deg from 119.5 W to 100.5 W:\n   ', ' '.join(f'{o:+.0f}/{h:.0f}' for o, h in zip(off, CORR_HALF)))
    print(f"  re-shaped from {fmt(n['overlay_from'])} ({n['overlay_km']:.0f} km to go, route/simulated x{n['overlay_stretch']:.4f}); the un-routed simulation ends at {n['sim_end'][1]:.3f} N {-n['sim_end'][0]:.3f} W")
    print(f"nominal vs official: transonic {fmt(n['transonic'])} vs {fmt(T['transonic'])}, subsonic {fmt(n['subsonic'])} vs {fmt(T['subsonic'])}, "
          f"landing burn {fmt(n['landing_burn'])} vs {fmt(T['lb'])}, catch {fmt(n['t_splash'])} vs {fmt(T['catch'])}")
