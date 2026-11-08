"""SGP4 TLE fitted to the Flight 15 model orbit (insertion to the deorbit burn) for the window-open T0; the page re-times the epoch
and RAAN for any T0 in the window (same Earth-fixed track). Placeholder catalog number 99995, designator 26999B. Model, not observed."""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/sda/starship-ITF-14/work'); sys.path.insert(0, '/Users/mickey/launch/Starship-FLT15/work')
from sgp4.api import Satrec, WGS72, jday
from sgp4 import exporter
from scipy.optimize import least_squares
from skyfield.api import load
from skyfield.sgp4lib import TEME_to_ITRF
import ift14_fit as F
import flt15_orbit as O
ts = load.timescale(builtin=True); pl, T = O.orbit, O.T

def dt_parts(t):
    d = np.datetime64(t).astype('datetime64[s]').astype(object); return d.year, d.month, d.day, d.hour, d.minute, d.second
def model_eci(t0, met_min):
    gmst0 = np.radians(ts.utc(*dt_parts(t0)).gmst*15.0)
    u = pl.u_nominal(met_min); Om = pl.L0 + gmst0 + pl.Odot*met_min
    x = O.RA*(np.cos(Om)*np.cos(u) - np.sin(Om)*np.sin(u)*pl.ci); y = O.RA*(np.sin(Om)*np.cos(u) + np.cos(Om)*np.sin(u)*pl.ci)
    return np.c_[x, y, O.RA*np.sin(u)*pl.si]
def make_sat(p, jd, fr, satnum):
    n, inc, raan, u0 = p; s = Satrec()
    s.sgp4init(WGS72, 'i', satnum, jd + fr - 2433281.5, 0.0, 0.0, 0.0, 1e-5, 0.0, np.radians(inc), np.radians(u0 % 360), n*2*np.pi/1440.0, np.radians(raan % 360))
    return s
def cks(l): return str(sum(int(c) if c.isdigit() else (1 if c == '-' else 0) for c in l[:68]) % 10)

def fit_tle(date, hms, satnum=99995):
    t0 = np.datetime64(f'{date}T{hms}'); te = t0 + np.timedelta64(int(round(T['ins0'])), 's'); jd, fr = jday(*dt_parts(te))
    met = np.arange(T['ins0']/60, T['deorb0']/60, 2.0); r_ref = model_eci(t0, met); tmin = met - T['ins0']/60
    def resid(p):
        _, r, _ = make_sat(p, jd, fr, satnum).sgp4_array(np.full(len(tmin), jd), fr + tmin/1440.0); return (r - r_ref).ravel()
    gmst0 = ts.utc(*dt_parts(t0)).gmst*15
    p0 = [1440/pl.T_nodal, pl.inc, (np.degrees(pl.L0) + gmst0) % 360, np.degrees(pl.u_ins) % 360]
    sol = least_squares(resid, p0, x_scale=[1e-3, 0.1, 0.1, 0.1]); rms = np.sqrt((sol.fun**2).mean()*3)
    sat = make_sat(sol.x, jd, fr, satnum)
    l1, l2 = exporter.export_tle(sat); l1 = l1[:9] + f'{"26999B":<8s}' + l1[17:]; l1 = l1[:68] + cks(l1); l2 = l2[:68] + cks(l2)
    name = 'STARSHIP FLT-15 PRELAUNCH MODEL - NOT OBSERVED'
    mv = np.arange(T['ins0']/60, T['deorb0']/60, 1.0); _, r, _ = sat.sgp4_array(np.full(len(mv), jd), fr + (mv - T['ins0']/60)/1440.0)
    rr = np.array([TEME_to_ITRF(j, ri, np.zeros(3))[0] for j, ri in zip(jd + fr + (mv - T['ins0']/60)/1440.0, r)])
    lon = np.degrees(np.arctan2(rr[:, 1], rr[:, 0])); lat = np.degrees(np.arctan(rr[:, 2]/np.hypot(rr[:, 0], rr[:, 1])/(1 - F.EF2)))
    lonm, latm = pl.nominal(mv); sep = F.gc_dist_deg(lon, lat, lonm, latm)
    print(f'{name}: n={sol.x[0]:.5f} i={sol.x[1]:.3f} RAAN={sol.x[2] % 360:.3f} u_epoch={sol.x[3] % 360:.3f}  fit rms {rms:.1f} km; track vs model max {sep.max():.3f} deg; epoch {str(te)[11:19]}Z (insertion)')
    return dict(name=name, l1=l1, l2=l2, t0=f'{date}T{hms}Z', raan_rate_deg_per_min=360.98564736629/1440.0, fit_rms_km=round(float(rms), 1), track_max_sep_deg=round(float(sep.max()), 3))

if __name__ == '__main__':
    out = fit_tle(O.DATE, O.T0_OPEN)
    json.dump(out, open(O.ROOT/'work'/'tle.json', 'w'), indent=1)
    (O.ROOT/'flt15_tle.txt').write_text(f"{out['name']} T0 {O.DATE} {O.T0_OPEN}Z\n{out['l1']}\n{out['l2']}\n")
    print(out['l1']); print(out['l2'])
