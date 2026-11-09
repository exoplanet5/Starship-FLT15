"""Export the Flight 15 model to the launch-page template's mission data (PROFILE_FORMAT.md): one JSON per profile in export/,
the nominal return to Starbase and the four ocean contingencies, sharing the ascent, coast and orbit and the Super Heavy return.
Run after flt15_descents.py --solve and flt15_tle.py:  ~/.venvs/astro313/bin/python flt15_export.py"""
import sys, json, math, os, numpy as np
from datetime import datetime, timedelta, timezone
from skyfield.api import Loader, load_file
from skyfield.framelib import itrs
sys.path.insert(0, '/Users/mickey/sda/starship-ITF-14/work'); sys.path.insert(0, '/Users/mickey/launch/Starship-FLT15/work')
import ift14_fit as F, ift14_entry as E
import flt15_orbit as O, flt15_zones as Z, flt15_descents as D
orbit, T, fmt = O.orbit, O.T, O.fmt
C = D.build(); TR, SOL = C['tracks'], C['solutions']
TLE = json.load(open(O.ROOT/'work'/'tle.json'))
OUT = O.ROOT/'export'; OUT.mkdir(exist_ok=True)
PROFILE = dict(id='starship-275', label='Starship · 275 km · 30.6°', family='low-inclination')
R6 = lambda v: round(float(v), 6)

# ---------------- model tracks to WGS-84 ----------------
def model_to_wgs84(lon, lat_geod, alt_model):
    lat_c = np.degrees(np.arctan((1 - F.EF2)*np.tan(np.radians(lat_geod)))); lo2, la2, h2 = [], [], []
    for lo, lc, a in zip(lon, lat_c, alt_model):
        g = E.geodetic((F.RE + a)*np.array([np.cos(np.radians(lc))*np.cos(np.radians(lo)), np.cos(np.radians(lc))*np.sin(np.radians(lo)), np.sin(np.radians(lc))]))
        la2.append(g[0]); lo2.append(g[1]); h2.append(g[2])
    return np.array(lo2), np.array(la2), np.array(h2)
def grid(t_a, t_b, dt_max):
    n = max(1, int(math.ceil((t_b - t_a)/dt_max - 1e-9))); return t_a + np.arange(n + 1)*(t_b - t_a)/n, (t_b - t_a)/n
def model_track(t_a, t_b, dt_max, color, label, taper=False):
    t, dt = grid(t_a, t_b, dt_max); lon, lat = orbit.nominal(t/60); alt = orbit.alt_nominal(t/60)
    lo, la, h = model_to_wgs84(lon, lat, alt)
    if taper: h = h - h[0]*(1 - t/T['seco'])**2; h[0] = 0.02          # the ascent starts on the pad (20 m: the engines on the mount, so the first sample is not on the ellipsoid)
    return dict(met0=float(t_a), dt=float(dt), lon=[R6(v) for v in lo], lat=[R6(v) for v in la], alt=[round(max(float(v), 0.0), 4) for v in h], color=color, label=label)
def uniform_track(tr, dt_max, color, label, dashed=False):
    """resample a descent (cached, 2 s samples + the landing point) onto a uniform grid, linear in ECEF"""
    t = np.array(tr['met']); X = np.array([E.ecef(la, lo, h) for lo, la, h in zip(tr['lon'], tr['lat'], tr['alt'])])
    tt, dt = grid(t[0], t[-1], dt_max); lon, lat, alt = [], [], []
    for q in tt:
        i = int(np.clip(np.searchsorted(t, q) - 1, 0, len(t) - 2)); f = (q - t[i])/(t[i + 1] - t[i])
        la, lo, h = E.geodetic(X[i] + (X[i + 1] - X[i])*f); lon.append(R6(lo)); lat.append(R6(la)); alt.append(round(max(float(h), 0.0), 4))
    out = dict(met0=float(tt[0]), dt=float(dt), lon=lon, lat=lat, alt=alt, color=color, label=label)
    if dashed: out['dashed'] = True
    return out

ASCENT = model_track(0.0, T['seco'], 1.0, '#d95926', 'Ascent', taper=True)
COAST = model_track(T['seco'], T['ins0'], 5.0, '#f2b35a', 'Suborbital coast')
def orbit_track(t_end): return model_track(T['ins0'], t_end, 20.0, '#c77dff', 'Orbit 275 km')
# Super Heavy: the booster's own height is above the ground path; blend the ascent's WGS-84 offset at MECO down to zero at the catch
BP = O.solve_booster(); B = O.booster_track(BP, 1.0)
lo_b, la_b = orbit.path_lonlat(B['s']); off_meco = ASCENT['alt'][int(T['meco'])] - float(orbit.h_of_t(T['meco']))
alt_b = B['h'] + off_meco*(1 - (B['t'] - T['meco'])/(T['bcatch'] - T['meco']))
BOOSTER = dict(met0=float(B['t'][0]), dt=1.0, lon=[R6(v) for v in lo_b], lat=[R6(v) for v in la_b], alt=[round(max(float(v), 0.0), 4) for v in alt_b], color='#9fb3c8', label='Super Heavy return')
DESC_STYLE = {'nominal': ('#4fc3f7', 'Return to Starbase'), 'indian': ('#ff5252', 'Indian Ocean contingency'), 'hawaii1': ('#ffb84f', 'Early Hawaii contingency'),
              'chile': ('#7ad97a', 'Chile contingency'), 'hawaii2': ('#c792ea', 'Late Hawaii contingency (orbit 16)')}
DESC = {k: uniform_track(TR[k], 2.0, *DESC_STYLE[k]) for k in TR}
DASHED = {k: dict(DESC[k], dashed=True) for k in DESC}

# ---------------- Sun and Moon (DE440s, shared kernel) ----------------
def body_table(body, start, end, model, step=60):
    load = Loader(str(O.ROOT/'work'), verbose=False); ts = load.timescale(builtin=True); eph = load_file(os.path.expanduser('~/JPLDE/de440s.bsp'))
    n = int((end - start).total_seconds()//step) + 2; utc = [start + timedelta(seconds=step*k) for k in range(n)]
    vec = eph['earth'].at(ts.from_datetimes(utc)).observe(eph[body]).frame_xyz(itrs).km.T
    return dict(start_utc=start.strftime('%Y-%m-%dT%H:%M:%SZ'), dt=step, vectors=[[round(float(v), 1) for v in p] for p in vec], model=model)
T0_LO = datetime(2026, 11, 13, 22, 45, tzinfo=timezone.utc); T0_HI = datetime(2026, 11, 14, 0, 0, tzinfo=timezone.utc)
END_MAX = max(tr['t_splash'] for tr in TR.values())
SUN = body_table('sun', T0_LO - timedelta(minutes=11), T0_HI + timedelta(seconds=END_MAX + 120), 'DE440s geometric Sun centre; WGS84 occultation')
MOON = body_table('moon', T0_LO - timedelta(minutes=11), T0_HI + timedelta(seconds=END_MAX + 120), 'DE440s geometric Moon centre')

# ---------------- shared metadata ----------------
STACK = dict(total_m=123, burn_end=END_MAX + 1, segments=dict(booster=dict(len=71, w=9, kind='solid'), ship=dict(len=52, w=9, kind='starship')),
             states=[dict(until=T['hotstage'], segs=['booster', 'ship'], zoom=1, burn='liquid'), dict(segs=['ship'], zoom=1.9, burn='liquid')])
COMPANION = dict(track='booster', label='Super Heavy', powered_intervals=[[T['bb0'], T['bb1']], [T['blb'], T['bcatch']]],
                 stack=dict(total_m=71, burn_end=T['bcatch'], segments=dict(booster=dict(len=71, w=9, kind='solid')), states=[dict(segs=['booster'], zoom=1.0, burn='liquid')]))
CITIES = [['Houston', -95.37, 29.76, 1], ['San Antonio', -98.49, 29.42, 1], ['Dallas', -96.80, 32.78, 1], ['Austin', -97.74, 30.27, 1], ['New Orleans', -90.07, 29.95, 1],
          ['Monterrey', -100.31, 25.67, 1], ['Mexico City', -99.13, 19.43, 1], ['Guadalajara', -103.35, 20.67, 1], ['Los Angeles', -118.24, 34.05, 1], ['Phoenix', -112.07, 33.45, 1],
          ['Miami', -80.19, 25.76, 1], ['Havana', -82.38, 23.13, 1], ['Honolulu', -157.86, 21.31, 1], ['Perth', 115.86, -31.95, 1], ['Sydney', 151.21, -33.87, 1],
          ['Jakarta', 106.85, -6.21, 1], ['Manila', 120.98, 14.60, 1], ['Tokyo', 139.69, 35.69, 1], ['Santiago', -70.65, -33.45, 1], ['Lima', -77.03, -12.05, 1],
          ['Cape Town', 18.42, -33.93, 1], ['Mumbai', 72.88, 19.08, 1], ['Nairobi', 36.82, -1.29, 1], ['São Paulo', -46.63, -23.55, 1], ['Buenos Aires', -58.38, -34.60, 1],
          ['Johannesburg', 28.04, -26.20, 1], ['Lagos', 3.39, 6.45, 1], ['Cairo', 31.24, 30.04, 1], ['Riyadh', 46.72, 24.69, 1], ['Dhaka', 90.41, 23.81, 1],
          ['Brownsville', -97.50, 25.90, 2], ['McAllen', -98.23, 26.20, 2], ['Corpus Christi', -97.40, 27.80, 2], ['Laredo', -99.51, 27.51, 2], ['El Paso', -106.49, 31.76, 2],
          ['Ciudad Juárez', -106.42, 31.69, 3], ['Chihuahua', -106.08, 28.63, 2], ['Hermosillo', -110.96, 29.07, 2], ['Mazatlán', -106.41, 23.22, 2], ['La Paz', -110.31, 24.14, 2],
          ['Saltillo', -101.00, 25.42, 2], ['Monclova', -101.42, 26.90, 3], ['Nuevo Laredo', -99.52, 27.49, 3], ['Reynosa', -98.30, 26.09, 3], ['Matamoros', -97.50, 25.87, 3],
          ['Tampico', -97.85, 22.25, 2], ['Mérida', -89.62, 20.97, 2], ['Cancún', -86.85, 21.16, 2], ['Hilo', -155.09, 19.72, 2], ['Lihue', -159.37, 21.98, 2],
          ['Darwin', 130.84, -12.46, 2], ['Port Moresby', 147.18, -9.44, 2], ['Brisbane', 153.03, -27.47, 2], ['Auckland', 174.76, -36.85, 2], ['Papeete', -149.57, -17.54, 2],
          ['Antofagasta', -70.40, -23.65, 2], ['Saint-Denis', 55.45, -20.88, 2], ['Port Louis', 57.50, -20.16, 2], ['Durban', 31.02, -29.86, 2], ['Colombo', 79.86, 6.93, 2],
          ['Guam', 144.79, 13.47, 2], ['Majuro', 171.38, 7.09, 3], ['Honiara', 159.97, -9.43, 3], ['Suva', 178.44, -18.14, 3], ['Hobart', 147.33, -42.88, 3],
          ['South Padre Island', -97.17, 26.11, 3], ['Harlingen', -97.70, 26.19, 3], ['Port Isabel', -97.21, 26.07, 3], ['Kingsville', -97.86, 27.52, 3], ['Victoria', -97.00, 28.81, 3],
          ['Ciudad Victoria', -99.15, 23.74, 3], ['Piedras Negras', -100.52, 28.70, 3], ['Hermosillo', -110.96, 29.07, 2], ['Culiacán', -107.39, 24.80, 3], ['Torreón', -103.44, 25.54, 3]]
seen = set(); CITIES = [c for c in CITIES if not (c[0] in seen or seen.add(c[0]))]
VIEWS = [['Global', -110, 0, 3.6], ['Starbase', -97.155, 25.7, 1.10], ['Approach', -102, 27.5, 1.45], ['Pacific entry', -150, 20, 2.4], ['Deorbit burn', 60, -25, 2.6]]
NOTES_COMMON = [
    'Pre-launch model, not telemetry and not an official ephemeris. Window 22:45-00:00 UTC on 13-14 Nov 2026 (16:45-18:00 CST), from the SpaceX mission page; T0 22:45 is the window open, not an announced liftoff time.',
    "Orbit plane: Flight 14's fit to the September 2026 navigational warnings (i 30.58 deg, 275 km circular with J2 secular rates, node 2.1 deg east of the pad, the ascent steering into the plane). SpaceX: Flight 15 flies the same altitude and inclination. The orbital rate is scaled by %+.4f (nodal period %.3f min, as a %.0f km mean orbit) so that the modeled catch happens at the official T+22:32:34; Flight 14's deorbit telemetry read 270 km against that model's 275." % (orbit.period_scale, orbit.T_nodal, O.RA*(1/(1 + orbit.period_scale))**(2/3) - F.RE),
    'Ascent: a shaped kinematic profile (Max Q 11.5 km at T+58 s, MECO 66 km and 1.5 km/s at T+2:19, SECO 150 km and 15 deg downrange at T+8:09, the Flight 14 assumptions) that ends on the suborbital coast ellipse in position and velocity. Insertion at T+26:10 is modeled at the burn start as an instantaneous circularisation at apogee.',
    'Super Heavy: 2-D point mass from the MECO state; 36 s boostback (%.0f m/s^2, %.0f m/s), ballistic fall with US76 drag (ballistic coefficient 3400 kg/m^2), 25 s landing burn to rest on the Pad 2 tower, arriving from the Gulf. Apogee %.0f km, %.0f km maximum downrange, %.0f m/s at the landing-burn start. Illustrative: the boostback and landing-burn throttles are not public.' % (B['a_bb'], B['dv_bb'], B['apogee'][1], B['max_downrange'][1], B['v_lb']*1000),
    "Descents: the Starship calibrated on Flight 14's reconstructed entry (3-DOF point mass, rotating WGS-84 Earth, J2, US Standard Atmosphere 1976, modified-Newtonian belly-first aerodynamics with a 55 deg hypersonic angle of attack, L/D 0.70, ballistic coefficient 331 kg/m^2, bank steering that fades out below Mach 1.5, a 10 s one-Raptor deorbit burn). The constant bank is a model choice; SpaceX describes a series of small turns over sparsely populated areas.",
    'Hazard areas and windows are the Flight 15 navigational warnings issued 7-9 Nov 2026: NAVAREA IV 1080/26 (launch areas A and B, 2245-0040Z), NAVAREA IV 1081/26 (ship catch: Laguna Madre, Padre Island and the Gulf off Starbase, 2108-0014Z on 14-15 Nov), HYDROPAC 3244/26 (Indian Ocean, 2253-0114Z), HYDROPAC 3257/26 (North Pacific, Midway Islands: early Hawaii, 0130-0501Z), HYDROPAC 3254/26 (South Pacific west of Chile, 0812-1137Z), HYDROPAC 3256/26 (North Pacific, Micronesia to north-east of Hawaii: late Hawaii, 2217-0144Z). The 14-20 Nov windows belong to the 13-19 Nov launch dates. Each contingency window opens as the modeled descent enters its area at about 80 km and closes about one orbit plus 20 min after the first-pass splashdown, so each covers two consecutive passes; the model shows the first.',
    'Sunlight: geometric WGS-84 occultation of the Sun centre with JPL DE440s at 60 s cadence. Shifting T0 keeps the Earth-fixed path and moves the Sun. Cities beyond the template list are added for the Americas, the Pacific and the deorbit regions; the sky chart opens on South Padre Island.',
    'Sources: spacex.com/launches/starship-flight-15 and spacex.com/updates (7 Nov 2026) for the timeline, window, catch description and contingency areas; the Flight 14 model at github.com/exoplanet5/Starship-IFT14 for the plane, zones and ship calibration.']
def ev(met, label, major=False, track=None):
    e = dict(met=float(met), label=label, major=major, kind='ship')
    if track: e['track'] = track
    return e
BOOSTER_EV = {'bb0', 'bb1', 'blb', 'bcatch'}
MAJOR = {'liftoff', 'meco', 'bcatch', 'seco', 'ins0', 'dep0', 'deorb0', 'entry', 'catch'}
KEYOF = {v: k for k, v in O._KEY.items()}
def official_events(until):
    out = []
    for t, txt in O.OFFICIAL:
        k = next((kk for kk, nm in O._KEY.items() if txt.startswith(kk)), None); key = O._KEY.get(k)
        if t <= until + 1e-6: out.append(ev(t, txt.split(' (')[0], key in MAJOR, 'booster' if key in BOOSTER_EV else None))
    return out
def descent_events(tr, key, land_label, burn_label=None):
    out = []
    if tr['t_burn'] is not None and burn_label: out += [ev(tr['t_burn'], burn_label, True), ev(tr['t_burn'] + D.BURN_S, 'Contingency deorbit burn shutdown (model)')]
    out += [ev(tr['t_ei'], '120 km (model)'), ev(tr['peak_heating'], f"Peak heating, {tr['peak_heating_alt']:.0f} km (model)"), ev(tr['t_max_g'], f"Max deceleration {tr['max_g']:.1f} g (model)"),
            ev(tr['transonic'], 'Transonic, Mach 1.2 (model)'), ev(tr['subsonic'], 'Subsonic (model)'), ev(tr['landing_burn'], 'Landing burn start (model)'), ev(tr['t_splash'], land_label, True)]
    return out

def zones_for(hit, impact=None):
    out = []
    for z in Z.ZONES:
        q = dict(key=z['key'], id=z['id'], role=z['role'], label=z['label'], color=z['color'], t0=z['t0'], t1=z['t1'], ring=z['ring'], centroid=z['centroid'], desc=z['desc'])
        if z['key'] == hit and impact is not None: q['impact'] = dict(met=float(impact[0]), lon=R6(impact[1]), lat=R6(impact[2]))
        elif z['role'] == 'launch': q['impact_note'] = 'Ascent overflight'
        elif z['role'] == 'landing': q['impact_note'] = 'Nominal return'
        else: q['impact_note'] = 'Contingency only' if hit == 'catch' else 'Other profile'
        out.append(q)
    return out

def mission(mid, mission_label, menu_label, desc_key, tracks, segments, events, zones, powered, phases, orbit_end, about, notes, met_default=150):
    end = tracks[segments[-1]]['met0'] + (len(tracks[segments[-1]]['lon']) - 1)*tracks[segments[-1]]['dt']
    meta = dict(title=f'Starship Flight 15 · {mission_label}', vehicle='Starship Flight 15', mission=mission_label, date=O.DATE, pad=[O.STARBASE[0], O.STARBASE[1]], pad_alt=0.0, pad_name='Starbase Pad 2',
                t0_default=O.T0_OPEN, t0_min=O.T0_OPEN, t0_max=O.T0_CLOSE,
                t0_presets=[dict(label='Open 22:45', time='22:45:00'), dict(label='23:10', time='23:10:00'), dict(label='23:35', time='23:35:00'), dict(label='Close 00:00', time='23:59:59')],
                met_default=met_default, insertion_met=float(T['ins0']), period_s=float(orbit.T_nodal*60), orbit_revolutions=(orbit_end - T['ins0'])/(orbit.T_nodal*60),
                orbit_end_met=float(orbit_end), inc=round(float(orbit.inc), 4), perigee_km=O.H_ORBIT, apogee_km=O.H_ORBIT, az=100.0,
                s4_cutoff_met=float(T['seco']), cutoff_label='SECO', stack=STACK, companion=COMPANION, extra_cities=CITIES, observer_home='South Padre Island',
                local_time=dict(label='CST', offset_h=-6), model_badge='Pre-launch model', readme_url='https://github.com/exoplanet5/Starship-FLT15#limitations', track_window=dict(track='orbit', revolutions=2), zone_noun='hazard area', layer_labels=dict(falls='Booster · other descents', zones='Hazard areas'),
                flight_segments=segments, powered_intervals=powered, phases=phases, views=VIEWS,
                visibility=dict(extent=[-125, -75, 10, 40], grid_step_deg=0.5, time_step_s=5, levels=[10]), about=about, notes=notes)
    if globals().get('_NO_WINDOW'): meta.pop('track_window', None)
    data = dict(schema_version=1, meta=meta, tracks=tracks, zones=zones, events=sorted(events, key=lambda e: e['met']),
                tle_ref=dict(t0=TLE['t0'], raan_rate_deg_per_min=TLE['raan_rate_deg_per_min'], sets=[dict(name=TLE['name'], l1=TLE['l1'], l2=TLE['l2'])]), sun=SUN, moon=MOON)
    out = OUT/f'{mid}.json'; out.write_text(json.dumps(data, separators=(',', ':'), ensure_ascii=False))
    print(f'wrote {out.name} ({out.stat().st_size/1e3:.0f} kB): end {fmt(end)}, {len(events)} events, ' + ', '.join(f'{k} {len(v["lon"])}' for k, v in tracks.items()))
    return dict(id=mid, label=menu_label, file=f'data/{mid}.json')

POW_ASC = [[0.0, T['meco']], [T['hotstage'], T['seco']], [T['ins0'], T['ins1']]]
PH_ASC = [dict(until=0, label='Pre-launch'), dict(until=T['meco'], label='Super Heavy ascent'), dict(until=T['hotstage'], label='Hot staging'), dict(until=T['seco'], label='Starship ascent'),
          dict(until=T['ins0'], label='Suborbital coast'), dict(until=T['ins1'], label='Orbit insertion burn'), dict(until=T['dep0'], label='Orbit · checkout'), dict(until=T['dep1'], label='Starlink V3 deploy')]
def descent_phases(tr, final):
    return [dict(until=tr['t_burn'] + D.BURN_S, label='Deorbit burn'), dict(until=tr['t_ei'], label='Entry coast'), dict(until=tr['transonic'], label='Atmospheric entry'),
            dict(until=tr['landing_burn'], label='Approach · belly flop'), dict(until=tr['t_splash'] + 1, label=final)]
ABOUT = ["Pre-launch model: Flight 14's zone-fitted orbit with the Flight 15 timeline, a Super Heavy return and a 3-DOF entry to the Pad 2 tower; not telemetry.",
         'Hazard areas and windows are the Flight 15 navigational warnings of 7-9 Nov 2026 (NAVAREA IV 1080/26 and 1081/26, HYDROPAC 3244, 3254, 3256 and 3257/26).',
         'WGS84 geometry · JPL DE440s sunlight.']
registry = []
_NO_WINDOW = False
# ---- nominal: return to launch site
n = TR['nominal']
tracks = dict(ascent=ASCENT, coast=COAST, orbit=orbit_track(T['deorb0']), descent=DESC['nominal'], booster=BOOSTER, **{k: DASHED[k] for k in ['indian', 'hawaii1', 'chile', 'hawaii2']})
events = official_events(1e9) + [ev(n['t_ei'], '120 km (model)'), ev(n['peak_heating'], f"Peak heating, {n['peak_heating_alt']:.0f} km (model)"), ev(n['t_max_g'], f"Max deceleration {n['max_g']:.1f} g (model)"),
                                  ev(TR['indian']['t_splash'], 'Indian Ocean splashdown (contingency)', False, 'indian'), ev(TR['hawaii1']['t_splash'], 'Early Hawaii splashdown (contingency)', False, 'hawaii1'),
                                  ev(TR['chile']['t_splash'], 'Chile splashdown (contingency)', False, 'chile')]
powered = POW_ASC + [[T['deorb0'], T['deorb1']], [T['lb'], T['catch']]]
phases = PH_ASC + [dict(until=T['deorb0'], label='Orbit coast'), dict(until=T['deorb1'], label='Deorbit burn'), dict(until=n['t_ei'], label='Entry coast'), dict(until=T['transonic'], label='Atmospheric entry'),
                   dict(until=T['lb'], label='Approach · belly flop'), dict(until=T['catch'] + 1, label='Landing burn · catch')]
registry.append(mission('flt15-rtls', 'Ship catch at Starbase', 'Flight 15 · ship catch at Starbase', 'nominal', tracks, ['ascent', 'coast', 'orbit', 'descent'], events,
                        zones_for('catch', (T['catch'], O.TOWER[0], O.TOWER[1])), powered, phases, T['deorb0'], ABOUT, NOTES_COMMON + [
    "Return: deorbit burn %.1f m/s at T+21:35:57 over the southern Indian Ocean, bank %.1f deg toward the south (the orbital pass runs about 400 km north of Starbase). 120 km at %s, 80 km at the official 'entry' time, peak heating %.0f km, %.1f g maximum. The pass crosses the Kauai Channel at about 80 km and northern Mexico at 60-50 km." % (SOL['nominal']['dv']*1000, SOL['nominal']['sigma'], fmt(n['t_ei']), n['peak_heating_alt'], n['max_g']),
    "Final approach from %s (the last %.0f km): the ground path is re-shaped to SpaceX's description and the ship-catch hazard area NAVAREA IV 1081/26, whose northern tip is at 26.77 N on the Laguna Madre: over the lagoon 65 km north of the South Padre Island town heading south-east, out over the Gulf to 36 km north-east of the tower, a final turn south-west to the tower. The simulation's arc length is kept (horizontal spacing stretched x%.2f), so heights, speeds and times are unchanged; the simulation's own constant-bank path would reach the tower directly from the west-north-west over Harlingen." % (fmt(n['overlay_from']), D.START_KM, n['overlay_stretch'])]))
# ---- contingencies
def contingency(mid, mission_label, menu_label, key, burn_label, land_label, land_point, extra_notes):
    tr = TR[key]; t_burn = tr['t_burn'] if tr['t_burn'] is not None else T['ins0']
    if key == 'indian':
        segs = ['ascent', 'descent']; tracks = dict(ascent=ASCENT, descent=DESC['indian'], booster=BOOSTER, rtls=DASHED['nominal'])   # the descent starts at SECO
        events = official_events(T['seco']) + [ev(T['ins0'], 'Insertion burn skipped (contingency)', True)] + descent_events(tr, key, land_label)
        powered = POW_ASC[:2]; phases = PH_ASC[:5] + [dict(until=tr['t_ei'], label='Ballistic coast, no insertion'), dict(until=tr['transonic'], label='Atmospheric entry'), dict(until=tr['landing_burn'], label='Belly flop'), dict(until=tr['t_splash'] + 1, label='Splashdown')]
        ins, orbit_end_eff, revs = float(T['seco']), None, None
    else:
        segs = ['ascent', 'coast', 'orbit', 'descent']; tracks = dict(ascent=ASCENT, coast=COAST, orbit=orbit_track(t_burn), descent=DESC[key], booster=BOOSTER, rtls=DASHED['nominal'])
        events = official_events(t_burn) + descent_events(tr, key, land_label, burn_label)
        powered = POW_ASC + [[t_burn, t_burn + D.BURN_S], [tr['landing_burn'], tr['t_splash']]]
        phases = PH_ASC + [dict(until=t_burn, label='Orbit coast')] + descent_phases(tr, 'Splashdown')
        if t_burn < T['dep1']: phases = [p for p in phases if p['until'] <= t_burn or p['label'] not in ('Orbit · checkout', 'Starlink V3 deploy')]
        orbit_end_eff = t_burn
    zones = zones_for(Z.ZONES[[z['key'] for z in Z.ZONES].index({'indian': 'indian', 'hawaii1': 'npac1', 'chile': 'chile', 'hawaii2': 'late'}[key])]['key'], (tr['t_splash'], tr['splash'][0], tr['splash'][1]))
    about = [f'Contingency profile: {mission_label}. ' + ABOUT[0]] + ABOUT[1:]
    global _NO_WINDOW; _NO_WINDOW = key == 'indian'
    m = mission(mid, mission_label, menu_label, key, tracks, segs, events, zones, powered, phases, orbit_end_eff if orbit_end_eff else tracks['descent']['met0'] + (len(tracks['descent']['lon']) - 1)*tracks['descent']['dt'], about, NOTES_COMMON + extra_notes)
    return m
# the Indian case has no orbit: insertion_met is set to SECO by overriding after the fact is awkward, so handle it inside mission(): use orbit_end = flight end
s_i, s_h1, s_c, s_h2 = SOL['indian'], SOL['hawaii1'], SOL['chile'], SOL['hawaii2']
registry.append(contingency('flt15-indian', 'Contingency · Indian Ocean', 'Contingency · no insertion, Indian Ocean', 'indian', None, 'Indian Ocean splashdown (model)', None,
    ["No insertion burn: from SECO the ship coasts on the suborbital ellipse (perigee %.0f km) with the 26 Starlink V3 still aboard (+50 t assumed, ballistic coefficient %.0f kg/m^2) and steers with a %.0f deg bank, reversed at %s, to the middle of the Indian Ocean area HYDROPAC 3244/26 (the Flight 14 polygon with a slightly longer western end; its window opens at SECO); %.1f g maximum." % (O.RP_COAST - F.RE, s_i['beta'], s_i['sigma'], fmt(s_i['t_rev']), TR['indian']['max_g'])]))
registry.append(contingency('flt15-hawaii-early', 'Contingency · early Hawaii', 'Contingency · early return, N Pacific (orbit 2)', 'hawaii1', 'Contingency deorbit burn, orbit 2 (model)', 'North Pacific splashdown, early Hawaii (model)', None,
    ["Early Hawaii: Flight 14's as-flown profile re-timed by the 42 s later insertion: deorbit burn %.1f m/s at %s on orbit 2 (the Early Hawaii Return Check is at T+01:27:46), bank %.1f deg, splashdown at Flight 14's reported point 25.50 N 155.43 W, 450 km north of Oahu. HYDROPAC 3257/26 is the Flight 14 North Pacific polygon; its window opens at 0130Z, as this descent enters the area at 80 km after a 22:45Z launch, and also covers an orbit-3 descent." % (s_h1['dv']*1000, fmt(s_h1['t_burn']), s_h1['sigma'])]))
registry.append(contingency('flt15-chile', 'Contingency · W of Chile', 'Contingency · mid-mission, W of Chile (orbit 7)', 'chile', 'Contingency deorbit burn, orbit 7 (model)', 'South Pacific splashdown W of Chile (model)', None,
    ["Mid-mission Chile: Flight 14's planned profile re-timed: deorbit burn %.1f m/s at %s over Tibet (the Mid Chile Return Check is at T+07:15:35), lift vector up, landing %s after the burn as planned on Flight 14, at %.1f S %.1f W, %.0f km off the Chilean coast. HYDROPAC 3254/26 is wider than the Flight 14 polygon; its window opens at 0812Z, as this descent enters the area at 80 km after a 22:45Z launch, and also covers an orbit-8 descent." % (s_c['dv']*1000, fmt(s_c['t_burn']), fmt(TR['chile']['t_splash'] - s_c['t_burn'])[2:], -TR['chile']['splash'][1], -TR['chile']['splash'][0], F.gc_dist_deg(TR['chile']['splash'][0], TR['chile']['splash'][1], -71.6, -29.9)*111.2)]))
registry.append(contingency('flt15-hawaii-late', 'Contingency · late Hawaii', 'Contingency · late return, N Pacific (orbit 16)', 'hawaii2', 'Contingency deorbit burn, orbit 16 (model)', 'North Pacific splashdown NE of Hawaii, late (model)', None,
    ["Late Hawaii on orbit 16, the pass after the Starbase pass: HYDROPAC 3256/26 runs from Micronesia to north-east of Hawaii and its window opens one orbit after the nominal entry, so the Starbase catch is declined at the nominal deorbit time and the ship flies one more orbit. Deorbit burn %.1f m/s (the Flight 14 value) at %s, %.0f min after the nominal burn and %.1f min before the same point of the next orbit, bank %.1f deg, splashdown at Flight 14's reported point %.2f N %.2f W inside the area's eastern lobe, %.0f km north-east of Honolulu. The window also covers an orbit-17 descent. The Final Hawaii Return Check at T+20:41:55 precedes both the nominal burn and this one." % (s_h2['dv']*1000, fmt(s_h2['t_burn']), (s_h2['t_burn'] - T['deorb0'])/60, (T['deorb0'] + orbit.T_nodal*60 - s_h2['t_burn'])/60, s_h2['sigma'], TR['hawaii2']['splash'][1], -TR['hawaii2']['splash'][0], F.gc_dist_deg(TR['hawaii2']['splash'][0], TR['hawaii2']['splash'][1], -157.86, 21.31)*111.2)]))
json.dump(dict(profile=PROFILE, missions=registry), open(OUT/'registry.json', 'w'), indent=1)
print('registry:', [m['id'] for m in registry])
