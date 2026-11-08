"""Hazard zones for the Flight 15 page and maps: the Flight 14 polygons (../navwarning, read through ift14_fit) with their daily
windows re-timed from the Flight 14 launch window (12:15-13:30Z) to the Flight 15 window (22:45-00:00Z). SpaceX says the Flight 15
contingency areas are the precleared areas used on Flight 14; no Flight 15 warnings were published when this was written, so the
windows below are model values, not published ones, and the page labels them so."""
import sys, numpy as np
sys.path.insert(0, '/Users/mickey/sda/starship-ITF-14/work'); sys.path.insert(0, '/Users/mickey/launch/Starship-FLT15/work')
import ift14_fit as F
import flt15_orbit as O

T0_OPEN14, T0_CLOSE14 = 12*3600 + 15*60, 13*3600 + 30*60
T0_OPEN15, T0_CLOSE15 = 22*3600 + 45*60, 24*3600            # 22:45Z to 00:00Z (next day)
def _hm(s): return int(s[:2])*3600 + int(s[2:])*60
def retime(t0, t1, late=False):
    """Flight 14 window HHMM-HHMM -> Flight 15 HHMM-HHMM keeping the offsets from the window open and close.
    late=True keeps the same offsets but one orbit-count later is not needed: the orbit-15 pass uses the nominal deorbit timing."""
    a, b = _hm(t0), _hm(t1)
    if b < a: b += 86400
    a15 = (a - T0_OPEN14 + T0_OPEN15) % 86400; b15 = (b - T0_CLOSE14 + T0_CLOSE15) % 86400
    return f'{a15//3600:02d}{a15%3600//60:02d}', f'{b15//3600:02d}{b15%3600//60:02d}'

def hhmm(sec): return f'{sec//3600:02d}{sec%3600//60:02d}'

ZONES = []
def add(key, fkey, role, label, color, ident, desc, window):
    z = F.Z[fkey]; pts = [(float(x), float(y)) for x, y in z['pts']]
    if z['wrap']: pts = [(((x + 180) % 360) - 180, y) for x, y in pts]          # back to -180..180; the page re-wraps edge by edge
    if pts[0] != pts[-1]: pts.append(pts[0])
    lon = np.array([p[0] for p in pts]); lat = np.array([p[1] for p in pts])
    ul = np.where(lon < 0, lon + 360, lon) if z['wrap'] else lon
    cen = [float(((ul.mean() + 180) % 360) - 180), float(lat.mean())]
    ZONES.append(dict(key=key, id=ident, role=role, label=label, color=color, t0=window[0], t1=window[1], ring=[[round(x, 4), round(y, 4)] for x, y in pts],
                      centroid=[round(cen[0], 3), round(cen[1], 3)], desc=desc, fkey=fkey, published=False))
SUN, DANGER, WARN, OK, VIOLET = '#ffd54f', '#ff5252', '#ffb84f', '#7ad97a', '#c792ea'
add('launchA', 'launchA', 'launch', 'Launch · Gulf', SUN, 'NAVAREA IV 922/26 A (Flight 14)', 'Launch hazard area A, Gulf of America and Yucatán Channel; Flight 14 polygon', retime('1215', '1414'))
add('launchB', 'launchB', 'launch', 'Launch · area B', SUN, 'NAVAREA IV 922/26 B (Flight 14)', 'Launch hazard area B; Flight 14 polygon', retime('1215', '1414'))
add('indian', 'indian', 'contingency', 'Indian Ocean · no insertion', DANGER, 'HYDROPAC 2751/26 (Flight 14)', 'Suborbital contingency: no insertion burn, splashdown from the SECO coast; Flight 14 polygon', retime('1223', '1447'))
add('npac1', 'npac', 'contingency', 'N Pacific · early Hawaii', WARN, 'NAVAREA XII 657/26 (Flight 14)', 'Early-mission contingency: deorbit on orbit 2, splashdown north of Hawaii; Flight 14 polygon', retime('1427', '1833'))
add('chile', 'chile', 'contingency', 'S Pacific · W of Chile', OK, 'HYDROPAC 2750/26 (Flight 14)', 'Mid-mission contingency: deorbit on orbit 7, landing several hundred km west of Chile; Flight 14 polygon', retime('2107', '0108'))
import json as _json
_DESC = O.ROOT/'work'/'descents.json'
def _late_box():
    """model splashdown area for the late-Hawaii option: a box around the solved splash point, longer uprange (south-west);
    not a published warning"""
    c = _json.load(open(_DESC)) if _DESC.exists() else {}
    tr = c.get('tracks', {}).get('hawaii2'); sol = c.get('solutions', {}).get('hawaii2')
    if not tr: return None
    lon, lat = tr['splash']; pts = [(lon - 6.0, lat - 2.6), (lon + 3.0, lat - 0.9), (lon + 3.4, lat + 1.4), (lon - 5.6, lat + 1.2), (lon - 6.0, lat - 2.6)]
    t_open = int(T0_OPEN15 + sol['t_burn']) % 86400; t_close = int(T0_CLOSE15 + tr['t_splash'] + 60*60) % 86400
    return pts, (hhmm(t_open), hhmm(t_close))
_lb = _late_box()
if _lb:
    pts, win = _lb; lon = np.array([p[0] for p in pts]); lat = np.array([p[1] for p in pts])
    ZONES.append(dict(key='late', id='Model area (not a warning)', role='contingency', label='SW of Hawaii · late (model)', color=VIOLET, t0=win[0], t1=win[1],
                      ring=[[round(x, 4), round(y, 4)] for x, y in pts], centroid=[round(float(lon[:-1].mean()), 3), round(float(lat[:-1].mean()), 3)],
                      desc='Late-mission contingency on the 14th orbit: deorbit 15-30 min before the nominal burn on the same Pacific pass, splashdown south-west of Hawaii short of the islands; a model area, no warning has been published', fkey=None, published=False))

if __name__ == '__main__':
    for z in ZONES: print(f"{z['key']:8s} {z['label']:28s} {z['t0']}-{z['t1']}Z  centroid {z['centroid']}  {len(z['ring'])} pts  {z['id']}")
