"""Starship Flight 15 maps (SatObserver-MX / Blue Marble style, ift14_style): the global mission profile and the catch approach.
Export: 300 dpi, 3960 px on the long side."""
import sys, numpy as np
sys.path.insert(0, '/Users/mickey/sda/starship-ITF-14/work'); sys.path.insert(0, '/Users/mickey/launch/Starship-FLT15/work')
import ift14_fit as F
from ift14_style import *
import ift14_style as S
from matplotlib.patches import Polygon as MPoly
from matplotlib.lines import Line2D
import flt15_orbit as O, flt15_zones as Z, flt15_descents as D
orbit, T, fmt = O.orbit, O.T, O.fmt
C = D.build(); TR, SOL = C['tracks'], C['solutions']
T0 = np.datetime64('2026-11-13T22:45:00')
FS = 0.6                     # font scale: the style module was tuned for 22 in at 180 dpi; this exports 13.2 in at 300 dpi
VIOLET, INK, ORANGE, STEEL = '#c792ea', '#e8eaed', '#ffa64d', '#b0c4de'

def new_fig(W, H, title):
    fig = plt.figure(figsize=(W, H), facecolor=BG)
    fig.text(0.012, 0.978, title, color='#e8eaed', fontsize=12.5*FS*1.25, fontweight='bold', va='center', ha='left')
    ax = fig.add_axes([0.012, 0.045, 0.976, 0.915]); ax.set_facecolor('#05070a')
    for s in ax.spines.values(): s.set_edgecolor('#2a323c')
    return fig, ax
def lab(ax, x, y, text, color='#ffffff', size=9.5, **kw): return label(ax, x, y, text, color=color, size=size*FS, **kw)
def sqm(ax, lon, lat, col, text=None, dx=1.6, dy=0, ha='left', va='center', hollow=False, size=8):
    ax.plot(lon, lat, 's', ms=size*FS*1.1, mfc='none' if hollow else col, mec=col if hollow else (8/255, 10/255, 14/255, 0.9), mew=1.4 if hollow else 1, zorder=7)
    if text: lab(ax, lon + dx, lat + dy, text, ha=ha, va=va)
def starm(ax, lon, lat, col, ms=16, hollow=False): star(ax, lon, lat, col, ms=ms*FS*1.15, hollow=hollow)
def credit_outside(fig): fig.text(0.988, 0.012, 'Illustration: Mickey · base map NASA Blue Marble', color=FGDIM, fontsize=7.5, ha='right', va='bottom')
def sun_lines(ax, glon, glat, alt):
    ax.contour(glon, glat, alt, levels=[0], colors=[SUN], linewidths=1.1, linestyles='solid', zorder=3)
    ax.contour(glon, glat, alt, levels=[-12], colors=[ACCENT], linewidths=1.0, linestyles='solid', zorder=3)
def legend2(ax, items, loc='lower left', anchor=(0.0, 0.0)):
    handles = [Line2D([], [], color=c, lw=lw, ls=ls, marker=m, ms=7*FS*1.2, mfc=c if m != '*' or hollow is False else 'none', mec=c, linestyle=ls) for (t, c, lw, ls, m, hollow) in items]
    lg = ax.legend(handles, [i[0] for i in items], loc=loc, bbox_to_anchor=anchor, ncol=int(np.ceil(len(items)/2)), fontsize=8*FS*1.15, frameon=True, facecolor=(0.04, 0.05, 0.07, 0.82),
                   edgecolor='#2a323c', labelcolor=FG, handlelength=2.2, columnspacing=1.1, borderpad=0.6)
    lg.set_zorder(12); return lg

def nominal_arrays():
    mm = np.arange(0, T['deorb0']/60 + 1e-9, 1/6); lo, la = orbit.nominal(mm); return mm*60, lo, la

# ================= global mission map =================
def unwrap(lon, center):
    """longitudes into [center-180, center+180) for a map centred on `center`"""
    lon = np.asarray(lon, float); return ((lon - center + 180) % 360) - 180 + center
def draw_base_centered(ax, t, center, grat=30):
    """draw_base for a world map centred on any meridian: the Blue Marble columns are rolled so the left edge is center-180"""
    img = S.earth_img(); Ww, Hh = img.size
    arr = np.asarray(img.reduce(3)); left = ((center - 180 + 180) % 360) - 180; c0 = int(round((left + 180)/360*arr.shape[1]))
    arr = np.roll(arr, -c0, axis=1); x0, x1 = center - 180, center + 180
    ax.imshow(arr, extent=[x0, x1, -90, 90], aspect='auto', interpolation='lanczos', zorder=0)
    ax.set_xlim(x0, x1); ax.set_ylim(-90, 90); ax.set_xticks([]); ax.set_yticks([])
    glon, glat = np.meshgrid(np.linspace(x0, x1, 1201), np.linspace(-90, 90, 801)); alt = S.sun_alt(t, glon, glat)
    sh = np.zeros(glon.shape + (4,)); sh[..., 2] = 10/255; sh[..., 3] = np.where(alt < 0, 0.35, 0)
    ax.imshow(sh, extent=[x0, x1, -90, 90], aspect='auto', origin='lower', zorder=1, interpolation='bilinear')
    slon, slat = S.sun_subpoint(t); slon = float(unwrap(slon, center))
    ax.plot(slon, slat, 'o', ms=7*FS*1.2, mfc=SUN, mec=(60/255, 40/255, 0, 0.6), mew=1, zorder=6)
    for k in range(8):
        a = k*np.pi/4; ax.plot([slon + np.cos(a)*1.9, slon + np.cos(a)*3.0], [slat + np.sin(a)*1.9, slat + np.sin(a)*3.0], color=SUN, lw=1.2, zorder=6)
    for lon in np.arange(np.ceil(x0/grat)*grat, x1 + 1e-9, grat): ax.axvline(lon, color=(1, 1, 1, 0.22 if lon % 360 == 0 else 0.10), lw=1, zorder=3)
    for lat in np.arange(-60, 61, grat): ax.axhline(lat, color=(1, 1, 1, 0.22 if lat == 0 else 0.10), lw=1, zorder=3)
    def lonlab(l):
        l = ((l + 180) % 360) - 180; return '0°' if l == 0 else ('180°' if abs(l) == 180 else f'{abs(int(l))}°{"W" if l < 0 else "E"}')
    def latlab(l): return '0°' if l == 0 else f'{abs(int(l))}°{"S" if l < 0 else "N"}'
    dx, dy = 360*0.0022, 180*0.0045
    for lon in np.arange(np.ceil(x0/grat)*grat, x1 - grat*0.5, grat): ax.text(lon + dx, -90 + dy, lonlab(lon), color=(1, 1, 1, 0.55), fontsize=9*FS*1.2, va='bottom', ha='left', path_effects=[HALO], zorder=8)
    for lat in np.arange(-60, 61, grat): ax.text(x0 + dx, lat + dy, latlab(lat), color=(1, 1, 1, 0.55), fontsize=9*FS*1.2, va='bottom', ha='left', path_effects=[HALO], zorder=8)
    return glon, glat, alt

def global_map(center=None):
    """the world map of the mission; center=None is the usual Greenwich-centred map, center=-150 puts the Pacific in the middle"""
    Wf = (lambda x: x) if center is None else (lambda x: unwrap(x, center))
    fig, ax = new_fig(13.2, 7.75, 'Starship Flight 15 mission profile · 2026-11-13 22:45 UTC (window 22:45–00:00Z) · ship catch at Starbase after ~14 orbits' + ('' if center is None else ' · Pacific view'))
    if center is None: glon, glat, alt = draw_base(ax, (-180, 180, -90, 90), T0)
    else: glon, glat, alt = draw_base_centered(ax, T0, center)
    sun_lines(ax, glon, glat, alt)
    for z in Z.ZONES:
        pts = np.array(z['ring']); lon = pts[:, 0].copy()
        if center is None:
            if lon.max() - lon.min() > 180: lon = np.where(lon < 0, lon + 360, lon)
            shifts = [0, -360] if lon.max() > 180 else [0]
        else: lon = Wf(lon); shifts = [0]
        for shift in shifts:
            q = np.c_[lon + shift, pts[:, 1]]
            ax.add_patch(MPoly(q, closed=True, facecolor=z['color'], alpha=0.20, edgecolor='none', zorder=4)); ax.add_patch(MPoly(q, closed=True, facecolor='none', edgecolor=z['color'], alpha=0.95, lw=1.0, zorder=4))
    ZL = {'launchA': (-92, 30.2, 'left', 'bottom', f"Launch · Gulf  {Z.ZONES[0]['t0']}–{Z.ZONES[0]['t1']}Z"), 'indian': (60, -34.5, 'left', 'top', f"Indian Ocean · no insertion  {Z.ZONES[2]['t0']}–{Z.ZONES[2]['t1']}Z"),
          'npac1': (178, 34.5, 'right', 'bottom', f"North Pacific · early Hawaii  {Z.ZONES[3]['t0']}–{Z.ZONES[3]['t1']}Z"),
          'late': (-178, 9.5, 'left', 'top', f"SW of Hawaii · late, orbit 14  {Z.ZONES[5]['t0']}–{Z.ZONES[5]['t1']}Z (model area)"),
          'chile': (-128, -18.0, 'left', 'bottom', f"S Pacific W of Chile · mid-mission  {Z.ZONES[4]['t0']}–{Z.ZONES[4]['t1']}Z")}
    if center is not None: ZL['npac1'] = (-178, 34.5, 'center', 'bottom', ZL['npac1'][4]); ZL['indian'] = (72, -31.5, 'left', 'top', ZL['indian'][4])
    for k, (x, y, ha, va, txt) in ZL.items():
        z = next(q for q in Z.ZONES if q['key'] == k); lab(ax, float(Wf(x)), y, txt, color=z['color'], ha=ha, va=va, linespacing=1.25)
    t, lo, la = nominal_arrays(); lo = Wf(lo)
    a = t < T['seco']; c = (t >= T['seco']) & (t <= T['ins0']); o = t >= T['ins0']
    plot_track(ax, lo[a], la[a], color=ORANGE, lw=2.2, alpha=0.95, zorder=5); plot_track(ax, lo[c], la[c], color=ORANGE, lw=1.6, ls=(0, (3, 3)), alpha=0.95, zorder=5)
    plot_track(ax, lo[o], la[o], color=STEEL, lw=0.8, alpha=0.7, zorder=5)
    for k, col, lw, ls in [('indian', DANGER, 1.5, (0, (3, 3))), ('hawaii1', WARN, 1.5, (0, (3, 3))), ('chile', OK, 1.5, (0, (3, 3))), ('hawaii2', VIOLET, 1.5, (0, (3, 3)))]:
        tr = TR[k]; plot_track(ax, Wf(tr['lon']), tr['lat'], color=col, lw=lw, ls=ls, alpha=0.95, zorder=6)
    nm = TR['nominal']; mm = np.array(nm['met']); hi = mm <= nm['t_ei']
    plot_track(ax, Wf(np.array(nm['lon'])[hi]), np.array(nm['lat'])[hi], color=INK, lw=1.1, ls=(0, (3, 3)), alpha=0.9, zorder=6)
    plot_track(ax, Wf(np.array(nm['lon'])[~hi]), np.array(nm['lat'])[~hi], color=INK, lw=2.2, alpha=0.95, zorder=6)
    i = np.where((np.diff(np.sign(la)) > 0) & (t[:-1] > T['ins0']))[0]
    for k, j in enumerate(i, 1): ax.text(lo[j] - 1.5, la[j] + 1.2, f'{k}', color=STEEL, fontsize=8*FS*1.2, ha='center', va='bottom', path_effects=[HALO], zorder=9)
    SBx = float(Wf(F.STARBASE[0]))
    ax.plot(SBx, F.STARBASE[1], marker='^', ms=10*FS*1.2, mfc=DANGER, mec=(8/255, 10/255, 14/255, 0.9), mew=1, ls='none', zorder=7)
    lab(ax, SBx + 1.0, F.STARBASE[1] + 2.0, 'Starbase · liftoff T+00:00:00 · booster catch T+00:06:35\nship catch T+22:32:34 on the orbit-15 pass', ha='left', va='bottom')
    j = int(np.searchsorted(t, T['seco'])); sqm(ax, lo[j], la[j], ORANGE, f'SECO {fmt(T["seco"])}', dx=1.8, dy=-1.2, va='top', hollow=True)
    j = int(np.searchsorted(t, T['ins0'])); sqm(ax, lo[j], la[j], ORANGE, f'insertion burn {fmt(T["ins0"])}', dx=2.0, dy=-0.8, va='top')
    j = int(np.searchsorted(t, T['dep0'])); sqm(ax, lo[j], la[j], STEEL, f'Starlink V3 deploy\n{fmt(T["dep0"])}–{fmt(T["dep1"])[2:]}', dx=-2.2, dy=-1.5, ha='right', va='top', hollow=True)
    n = TR['nominal']
    sqm(ax, float(Wf(n['burn'][0])), n['burn'][1], INK, f'deorbit burn {fmt(T["deorb0"])}', dx=2.0, dy=1.5, va='bottom')
    sqm(ax, float(Wf(n['ei'][0])), n['ei'][1], INK, f'120 km {fmt(n["t_ei"])}', dx=-2.0, dy=1.4, ha='right', va='bottom', hollow=True)
    jj = int(np.argmin(np.abs(np.array(n['met']) - T['entry']))); sqm(ax, float(Wf(n['lon'][jj])), n['lat'][jj], INK, f'entry {fmt(T["entry"])}', dx=-1.5, dy=2.0, ha='right', va='bottom', hollow=True)
    for k, col, txt, dx, dy, ha, va in [('indian', DANGER, f"no insertion → splash {fmt(TR['indian']['t_splash'])}", 2.5, -1.5, 'left', 'top'), ('hawaii1', WARN, f"orbit 2 → splash {fmt(TR['hawaii1']['t_splash'])}", 2.5, -2.2, 'left', 'top'),
                                         ('chile', OK, f"orbit 7 → landing {fmt(TR['chile']['t_splash'])}", 2.5, 2.0, 'left', 'bottom'), ('hawaii2', VIOLET, f"orbit 14 → splash {fmt(TR['hawaii2']['t_splash'])}", 2.5, -2.2, 'left', 'top')]:
        tr = TR[k]; sx = float(Wf(tr['splash'][0])); starm(ax, sx, tr['splash'][1], col, ms=14, hollow=True); lab(ax, sx + dx, tr['splash'][1] + dy, txt, ha=ha, va=va)
        flip = center is not None and k in ('hawaii1', 'hawaii2')       # these burns sit near the Pacific view's right edge
        below = k == 'hawaii1' and not flip                               # the orbit-2 burn label goes under its square, except in the Pacific view
        sqm(ax, float(Wf(tr['burn'][0])), tr['burn'][1], col, f"burn {fmt(tr['t_burn'])}" if tr['t_burn'] else None, dx=(-1.8 if flip else 1.8), ha=('right' if flip else 'left'), dy=(-1.5 if below else 1.2), va=('top' if below else 'bottom'))
    starm(ax, SBx, F.STARBASE[1], INK, ms=18)
    hud(ax, [('STARSHIP FLIGHT 15 · first ship catch', FG, True),
             (f'i = {orbit.inc:.2f}°   h = 275 km   nodal period {orbit.T_nodal:.2f} min', FG, False),
             (f'launch 13 Nov 2026 22:45–00:00Z (T+ from 22:45Z) · Pad 2', FG, False),
             (f'MECO {fmt(T["meco"])}  booster catch {fmt(T["bcatch"])}  SECO {fmt(T["seco"])}  insertion {fmt(T["ins0"])}', SUN, False),
             (f'return: deorbit {fmt(T["deorb0"])} ({SOL["nominal"]["dv"]*1000:.0f} m/s)  entry {fmt(T["entry"])}  catch {fmt(T["catch"])}', INK, False),
             (f'Indian Ocean   no insertion burn        splash {fmt(TR["indian"]["t_splash"])}', DANGER, False),
             (f'N Pacific      deorbit {fmt(TR["hawaii1"]["t_burn"])} (orbit 2)  splash {fmt(TR["hawaii1"]["t_splash"])}', WARN, False),
             (f'W of Chile     deorbit {fmt(TR["chile"]["t_burn"])} (orbit 7)  landing {fmt(TR["chile"]["t_splash"])}', OK, False),
             (f'SW of Hawaii   deorbit {fmt(TR["hawaii2"]["t_burn"])} (orbit 14) splash {fmt(TR["hawaii2"]["t_splash"])}', VIOLET, False),
             ('zones: Flight 14 warnings re-timed; late-Hawaii area from the model; none published', FGDIM, False)], x0=(0.010 if center is None else 0.57), y0=0.062, width=0.40, lh=0.024, fontsize=11*FS*1.1)
    legend2(ax, [('ascent', ORANGE, 2.2, '-', None, False), ('coast', ORANGE, 1.6, (0, (3, 3)), None, False), ('orbit', STEEL, 0.9, '-', None, False), ('return to Starbase', INK, 2.0, '-', None, False),
                 ('contingency descent', WARN, 1.5, (0, (3, 3)), None, False), ('Sun 0°', SUN, 1.1, '-', None, False), ('Sun −12°', ACCENT, 1.0, '-', None, False), ('splashdown', INK, 0, 'none', '*', True)],
            loc=('lower right' if center is None else 'upper left'), anchor=((1.0, 0.0) if center is None else (0.0, 1.0)))
    credit_outside(fig)
    out = O.ROOT/('flt15_mission_map.png' if center is None else 'flt15_mission_map_pacific.png'); fig.savefig(out, dpi=300, facecolor=BG); print('wrote', out)

# ================= catch approach map =================
def approach_map():
    fig, ax = new_fig(13.2, 8.4, 'Starship Flight 15 return to Starbase · 2026-11-14 21:17 UTC (catch, T0 22:45Z) · entry over the Pacific, approach over northern Mexico and south Texas')
    t_catch = T0 + np.timedelta64(int(T['catch']), 's')
    ext = (-120, -92, 20, 34); glon, glat, alt = draw_base(ax, ext, t_catch, grat=5, label_step=5); sun_lines(ax, glon, glat, alt)
    n = TR['nominal']; met = np.array(n['met']); lon = np.array(n['lon']); lat = np.array(n['lat']); al = np.array(n['alt']); sp = np.array(n['speed'])
    t_p = np.arange(T['deorb0'], T['catch'] + 600, 10.0); lo_p, la_p = orbit.lonlat(t_p/60, orbit.u_nominal(t_p/60)); plot_track(ax, lo_p, la_p, color=STEEL, lw=0.9, alpha=0.7, zorder=5)   # the pass, had the ship stayed in orbit
    plot_track(ax, lon, lat, color=INK, lw=2.2, alpha=0.95, zorder=6)
    for tm in range(int(T['entry']) + 2*60, int(T['catch']) - 60, 120):
        j = int(np.argmin(np.abs(met - tm)))
        if ext[0] < lon[j] < ext[1] and ext[2] < lat[j] < ext[3]:
            ax.plot(lon[j], lat[j], 'o', ms=4*FS*1.2, mfc=INK, mec='#05070a', zorder=7)
            rho, Tk, c = D.atmo(al[j]); lab(ax, lon[j] + 0.25, lat[j] - 0.35, f'{fmt(tm)[2:]}  {al[j]:.0f} km  M{sp[j]*1000/c:.0f}', size=8, va='top')
    for name, (la_, lo_) in {**D.PLACES, 'San Antonio': (29.42, -98.49), 'Houston': (29.76, -95.37), 'Mazatlán': (23.22, -106.41), 'La Paz': (24.14, -110.31), 'Torreón': (25.54, -103.44), 'Culiacán': (24.80, -107.39), 'Saltillo': (25.42, -101.0)}.items():
        if name in ('Lihue (Kauai)', 'Honolulu', 'Hilo', 'South Padre Island'): continue
        if ext[0] < lo_ < ext[1] and ext[2] < la_ < ext[3]:
            ax.plot(lo_, la_, 'o', ms=3.2*FS*1.2, mfc='#cfd6dd', mec='#05070a', mew=0.6, zorder=7); lab(ax, lo_ + 0.18, la_ + 0.12, name, size=8, color='#e8eaed', va='bottom')
    ax.plot(*F.STARBASE, marker='^', ms=10*FS*1.2, mfc=DANGER, mec=(8/255, 10/255, 14/255, 0.9), mew=1, ls='none', zorder=8)
    lab(ax, F.STARBASE[0] - 0.35, F.STARBASE[1] - 0.45, f'Starbase · catch {fmt(T["catch"])}', ha='right', va='top')
    # inset: final approach
    ins = ax.inset_axes([0.015, 0.055, 0.37, 0.43]); ins.set_facecolor('#05070a'); ins.set_zorder(20)
    iext = (-97.70, -96.85, 25.72, 26.40); S.draw_base(ins, iext, t_catch, grat=0.5, night=False, sun=False, label_step=10)
    k = np.where((met >= n['overlay_from'] - 150))[0]; ins.plot(lon[k], lat[k], color=INK, lw=2.4, zorder=6); ins.set_xlim(iext[0], iext[1]); ins.set_ylim(iext[2], iext[3])
    for tm in [T['transonic'], T['subsonic'], T['lb']]:
        j = int(np.argmin(np.abs(met - tm))); ins.plot(lon[j], lat[j], 's', ms=5*FS*1.2, mfc='none', mec=INK, mew=1.2, zorder=7)
        ins.text(lon[j] + 0.012, lat[j] + 0.006, f'{fmt(tm)[2:]}  {al[j]:.0f} km', color='#ffffff', fontsize=7.5*FS*1.15, path_effects=[HALO], zorder=9)
    ins.plot(*F.STARBASE, marker='^', ms=9*FS*1.2, mfc=DANGER, mec='#05070a', zorder=8)
    for name, (la_, lo_) in {'South Padre Island': (26.11, -97.17), 'Port Isabel': (26.07, -97.21), 'Brownsville': (25.90, -97.50), 'Laguna Madre': (26.30, -97.35)}.items():
        ins.plot(lo_, la_, 'o', ms=3*FS*1.2, mfc='#cfd6dd', mec='#05070a', mew=0.5, zorder=7); ins.text(lo_ + 0.012, la_ + 0.008, name, color='#e8eaed', fontsize=7.5*FS*1.15, path_effects=[HALO], zorder=9)
    ins.text(0.02, 0.97, 'final approach: north of South Padre Island,\nout over the Gulf, final turn south to the tower', transform=ins.transAxes, color=FG, fontsize=7.5*FS*1.15, va='top', path_effects=[HALO], zorder=9)
    for s in ins.spines.values(): s.set_edgecolor('#8899aa')
    hud(ax, [('RETURN TO STARBASE · model', FG, True),
             (f'deorbit burn {fmt(T["deorb0"])}  {SOL["nominal"]["dv"]*1000:.0f} m/s, bank {SOL["nominal"]["sigma"]:.0f}° south', FG, False),
             (f'120 km {fmt(n["t_ei"])}   entry (80 km) {fmt(T["entry"])}   peak heating {fmt(n["peak_heating"])} {n["peak_heating_alt"]:.0f} km', FG, False),
             (f'max {n["max_g"]:.1f} g   transonic {fmt(T["transonic"])}   landing burn {fmt(T["lb"])}   catch {fmt(T["catch"])}', FG, False),
             ('orbital pass ~400 km N of Starbase; the glide banks south', FGDIM, False),
             (f'last {D.START_KM:.0f} km re-shaped to the described approach', FGDIM, False)], x0=0.012, y0=0.495, width=0.37, lh=0.024, fontsize=11*FS*1.1)
    legend2(ax, [('return glide (2-min marks: MET, height, Mach)', INK, 2.2, '-', None, False), ('orbital pass, no deorbit', STEEL, 0.9, '-', None, False), ('Sun 0°', SUN, 1.1, '-', None, False), ('Sun −12°', ACCENT, 1.0, '-', None, False)], loc='upper left', anchor=(0.0, 1.0))
    credit_outside(fig)
    out = O.ROOT/'flt15_catch_approach.png'; fig.savefig(out, dpi=300, facecolor=BG); print('wrote', out)

if __name__ == '__main__':
    global_map(); global_map(center=-150); approach_map()
