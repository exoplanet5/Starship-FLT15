"""Hazard areas for the Flight 15 page and maps: the Flight 15 navigational warnings in ../navwarning (NGA text as issued
7-9 Nov 2026) parsed into polygons, with their published daily windows (HHMM UTC).

  NAVAREA IV 1080/26 A, B   launch, Gulf of America and Yucatan Channel; booster return      2245-0040Z 13 Nov, daily 14-19 Nov
  NAVAREA IV 1081/26        ship catch: Laguna Madre, Padre Island and the Gulf off Starbase  2108-0014Z 14-15 Nov, daily 15-20 Nov
  HYDROPAC 3244/26          eastern Indian Ocean, no-insertion contingency                    2253-0114Z daily 13-19 Nov
  HYDROPAC 3257/26          North Pacific, Midway Islands: early Hawaii, orbits 2-3           0130-0501Z daily 14-20 Nov (the Flight 14 polygon)
  HYDROPAC 3254/26          South Pacific west of Chile: mid-mission, orbits 7-8              0812-1137Z daily 14-20 Nov
  HYDROPAC 3256/26          North Pacific, Micronesia to north-east of Hawaii: late Hawaii,   2217-0144Z daily 14-20 Nov
                            the pass after the Starbase pass (orbit 16) and the next one
The dates are the UTC days of each window's opening; the launch window is 13 Nov 22:45-00:00Z, so the 14-20 Nov windows
belong to launches on 13-19 Nov. Each contingency window opens as the modeled descent enters the area and closes about
one orbit plus 20 minutes after the first-pass splashdown, i.e. each covers two consecutive passes."""
import re, numpy as np
from pathlib import Path
from shapely.geometry import Polygon, LineString, Point

ROOT = Path('/Users/mickey/launch/Starship-FLT15'); NAV = ROOT/'navwarning'
T0_OPEN15, T0_CLOSE15 = 22*3600 + 45*60, 24*3600            # launch window 22:45Z to 00:00Z (next day)

def _dm(tok):
    """'24-58.00S' -> -24.9667, '070-56.00E' -> 70.9333"""
    d, mi, h = re.fullmatch(r'(\d{2,3})-(\d{2}(?:\.\d+)?)([NSEW])', tok).groups()
    v = int(d) + float(mi)/60; return -v if h in 'SW' else v
_PAIR = re.compile(r'(\d{2}-\d{2}\.\d{2}[NS])\s+(\d{3}-\d{2}\.\d{2}[EW])')

def read_warning(path):
    """one NGA text warning -> list of areas: dict(name, poly, pts, wrap, t0, t1, dates, issued, cancel)"""
    txt = Path(path).read_text(); lines = txt.splitlines()
    issued = lines[0].strip(); name = lines[1].strip().rstrip('.')
    m = re.search(r'(\d{2})(\d{4})Z TO (\d{2})(\d{4})Z NOV', txt)                      # dated form: 132245Z TO 140040Z NOV
    md = re.search(r'(\d{4})Z TO (\d{4})Z DAILY (\d+) THRU? (\d+) NOV', txt)           # daily form ('THR' in 3256/26)
    if m: t0, t1, dates = m.group(2), m.group(4), f'{int(m.group(1))} Nov'
    else: t0, t1, dates = md.group(1), md.group(2), ''
    if md: dates += (', ' if dates else '') + f'daily {int(md.group(3))}-{int(md.group(4))} Nov'
    cancel = re.search(r'CANCEL THIS MSG (\d{6})Z', txt).group(1)
    body = txt.split('BOUND BY', 1)[1]
    parts = re.split(r'\n\s+([A-Z])\.\s+', '\n' + body) if 'AREAS BOUND BY' in txt else [None, None, body]
    out = []
    for sub, chunk in zip(parts[1::2], parts[2::2]):
        pts = [(_dm(lo), _dm(la)) for la, lo in _PAIR.findall(chunk)]
        lons = np.array([p[0] for p in pts]); wrap = lons.max() - lons.min() > 180
        if wrap: pts = [(x % 360, y) for x, y in pts]                                   # keep the ring continuous in 0..360
        if pts[0] != pts[-1]: pts.append(pts[0])
        out.append(dict(name=name + (f' {sub}' if sub else ''), sub=sub, poly=Polygon(pts), pts=pts, wrap=wrap, t0=t0, t1=t1, dates=dates,
                        issued=issued, cancel=cancel, file=Path(path).name))
    return out

AREAS = {a['name']: a for f in sorted(NAV.iterdir()) if f.is_file() and not f.suffix and f.name != 'README.md' for a in read_warning(f)}
KEYS = {'launchA': 'NAVAREA IV 1080/26 A', 'launchB': 'NAVAREA IV 1080/26 B', 'indian': 'HYDROPAC 3244/26', 'npac1': 'HYDROPAC 3257/26',
        'chile': 'HYDROPAC 3254/26', 'late': 'HYDROPAC 3256/26', 'catch': 'NAVAREA IV 1081/26'}
ZD = {k: AREAS[v] for k, v in KEYS.items()}                  # ift14_fit-compatible zone dicts (poly, pts, wrap): F.zcontains, F.crossings work on them

def contains(z, lon, lat):
    p = Point(lon, lat); return bool(z['poly'].contains(p) or (z['wrap'] and z['poly'].contains(Point(lon % 360, lat))))
def crossings(z, met, lon, lat):
    """[(met_in, met_out), ...] intervals of a sampled track inside the area"""
    inside = np.array([contains(z, x, y) for x, y in zip(lon, lat)]); idx = np.where(inside)[0]
    return [] if len(idx) == 0 else [(float(met[g[0]]), float(met[g[-1]])) for g in np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)]
def centerline(z, step=0.5):
    """(lon, lat) mid-points of north-south cuts through the area, west to east"""
    x0, _, x1, _ = z['poly'].bounds; pts = []
    for lon in np.arange(x0 + step, x1, step):
        cut = z['poly'].intersection(LineString([(lon, -89), (lon, 89)]))
        if not cut.is_empty: pts.append((((lon + 180) % 360) - 180, 0.5*(cut.bounds[1] + cut.bounds[3])))
    return np.array(pts)
def hhmm(sec): return f'{sec//3600:02d}{sec%3600//60:02d}'
def window_met(z):
    """the published window as MET offsets (s) for the earliest and the latest launch: (open, close) for T0 22:45Z, (open, close) for T0 00:00Z;
    the window's first date gives its day relative to the 13 Nov launch date"""
    a = int(z['t0'][:2])*3600 + int(z['t0'][2:])*60; b = int(z['t1'][:2])*3600 + int(z['t1'][2:])*60
    if b < a: b += 86400
    day = int(re.search(r'(\d+)', z['dates']).group(1)) - 13; a += day*86400; b += day*86400
    return (a - T0_OPEN15, b - T0_OPEN15), (a - T0_CLOSE15, b - T0_CLOSE15)

SUN, DANGER, WARN, OK, VIOLET, CYAN = '#ffd54f', '#ff5252', '#ffb84f', '#7ad97a', '#c792ea', '#4fc3f7'
ZONES = []
def add(key, role, label, color, desc):
    z = ZD[key]; pts = [(float(x), float(y)) for x, y in z['pts']]
    if z['wrap']: pts = [(((x + 180) % 360) - 180, y) for x, y in pts]          # back to -180..180; the page re-wraps edge by edge
    lon = np.array([p[0] for p in pts]); lat = np.array([p[1] for p in pts])
    ul = np.where(lon < 0, lon + 360, lon) if z['wrap'] else lon
    cen = [float(((ul[:-1].mean() + 180) % 360) - 180), float(lat[:-1].mean())]
    ZONES.append(dict(key=key, id=z['name'], role=role, label=label, color=color, t0=z['t0'], t1=z['t1'], dates=z['dates'], issued=z['issued'],
                      ring=[[round(x, 4), round(y, 4)] for x, y in pts], centroid=[round(cen[0], 3), round(cen[1], 3)], desc=desc, published=True))
add('launchA', 'launch', 'Launch · Gulf', SUN, 'Launch hazard area A, Gulf of America, Yucatán Channel and the Caribbean south of Cuba (NAVAREA IV 1080/26, 9 Nov 2026)')
add('launchB', 'launch', 'Launch · area B', SUN, 'Launch hazard area B off Starbase, the Super Heavy return (NAVAREA IV 1080/26, 9 Nov 2026)')
add('indian', 'contingency', 'Indian Ocean · no insertion', DANGER, 'Suborbital contingency: no insertion burn, splashdown from the SECO coast (HYDROPAC 3244/26, 7 Nov 2026)')
add('npac1', 'contingency', 'N Pacific · early Hawaii', WARN, 'Early-mission contingency: deorbit on orbit 2 (the window also covers orbit 3), splashdown north of Hawaii; the Flight 14 polygon (HYDROPAC 3257/26, 9 Nov 2026)')
add('chile', 'contingency', 'S Pacific · W of Chile', OK, 'Mid-mission contingency: deorbit on orbit 7 (the window also covers orbit 8), landing several hundred km west of Chile (HYDROPAC 3254/26, 9 Nov 2026)')
add('late', 'contingency', 'N Pacific · late Hawaii', VIOLET, 'Late-mission contingency on orbit 16, the pass after the Starbase pass (the window also covers orbit 17): entry from Micronesia to a splashdown north-east of Hawaii (HYDROPAC 3256/26, 9 Nov 2026)')
add('catch', 'landing', 'Ship catch · Starbase', CYAN, 'Ship return: Laguna Madre, Padre Island and the Gulf off Starbase for the final approach and catch; the land corridor is under FAA notices (NAVAREA IV 1081/26, 9 Nov 2026)')

if __name__ == '__main__':
    import sys
    def fmt(s): s = int(round(s)); return f'T+{s//3600:02d}:{s%3600//60:02d}:{s%60:02d}'
    for z in ZONES:
        (a, b), (c, d) = window_met(ZD[z['key']])
        print(f"{z['key']:8s} {z['id']:22s} {z['t0']}-{z['t1']}Z {z['dates']:24s} early launch {fmt(a)}..{fmt(b)}  late launch {fmt(c)}..{fmt(d)}  {len(z['ring'])} pts  centroid {z['centroid']}")
