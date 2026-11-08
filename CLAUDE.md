# CLAUDE.md

Starship Flight 15 (NET 13 Nov 2026, window 22:45–00:00 UTC, Pad 2): the first ship catch at Starbase. Pre-launch model on
Flight 14's zone-fitted orbit and calibrated ship (`~/sda/starship-ITF-14`, imported, not copied), the official Flight 15
timeline, a Super Heavy return and a 3D page built from `~/launch/launch-page-template` with the `launch-page` skill.
Repo: github.com/exoplanet5/Starship-FLT15; GitHub Pages serves `main` → `/docs`, so a push to `main` deploys https://exoplanet5.github.io/Starship-FLT15/ (publish with the skill's `publish_page.py`, only when asked). README.md has the results, assumptions and sources. Everything runs with `~/.venvs/astro313/bin/python`.

## Pipeline (order matters)

1. `work/flt15_orbit.py`: plane (from `ift14_fit`), Flight 15 timeline (`flight-timeline.txt`), shaped ascent, booster. Reads `work/phase.json` for the period scale if it exists.
2. `work/flt15_descents.py --solve`: nominal (solves Δv, bank and the period scale → writes `phase.json`), then the four contingencies and the sampled tracks → `work/descents.json`. About 10 min. Without `--solve` it only fills in missing solutions and tracks.
3. `work/flt15_tle.py` → `work/tle.json`, `flt15_tle.txt` (after step 2: it uses the period scale).
4. `work/flt15_map.py` → the three PNGs (`global_map()`, `global_map(center=-150)` for the Pacific view, `approach_map()`). `work/flt15_zones.py` holds the zones and re-timed windows.
5. `work/flt15_export.py` → `export/*.json` and `export/registry.json`; then `new_page.py` (first mission without `--add`, the rest with `--add`); the menu order in `docs/data/profiles.json` was set by hand to nominal, Indian, early Hawaii, Chile, late Hawaii.

## Key facts and decisions

- **User inputs (8 Nov 2026):** same altitude and inclination as Flight 14; first full catch (booster and ship); more than 22 h to align with Starbase; contingencies as Flight 14 plus a late-mission Hawaii option; no Flight 15 warnings yet, adopt Flight 14's zones and re-time them; extend the model to reentry and catch.
- **Official timeline** from spacex.com/launches/starship-flight-15 (fetched 8 Nov via Exa; the page is client-rendered, plain fetches return an empty shell). Window 22:45–00:00Z = 16:45–18:00 CST. `t0_max` on the page is 23:59:59 because the template needs a same-date interval.
- **Period scale +0.00134** (nodal period 89.63 min, as a 269 km mean orbit) is the one free orbit parameter, solved with the nominal Δv and bank for the official catch time. Without it the catch came 1 min 46 s late while the altitude at the official entry time (80 km) already matched Flight 14's convention. Flight 14's deorbit telemetry read 270 km, so a lower mean orbit is the physical reading. Do not retune it by hand; re-solve the nominal.
- **Return pass:** the 15th pass after insertion crosses Starbase's longitude ~400 km north of the pad; the ship banks ~20° south, like Flight 14's reconstructed entry. The pass crosses the Kauai Channel at ~80 km (Honolulu 49 km, Lihue 81 km), Chihuahua 18 km at 61 km. Constant bank is a model choice.
- **Late Hawaii = 14th orbit, before the Starbase pass (user, 8 Nov).** On the same Pacific pass as the return: burn time and bank solved with a fixed 69 m/s for a splashdown at 18.3 N 165.5 W (short of the islands): T+21:19:46, 16 min before the nominal burn, lift up. The splashdown area is a model box built in `flt15_zones._late_box` from `descents.json`, not the Flight 14 North Pacific zone (user: need not follow the early-Hawaii zone). Rejected: one orbit later (user), and a shortened glide to the Flight 14 lobe beyond the islands (130° bank, 37 g).
- **Chile is solved to Flight 14's planned landing point** (29.9 S 86.9 W), not its landing time: with the period scale, a time-solved descent landed 3° east, outside the zone.
- **Early Hawaii** is Flight 14's as-flown burn re-timed by the 42 s later insertion (T+02:12:42) and solved to Flight 14's reported splashdown.
- **Indian Ocean** as on Flight 14: coast from SECO with +50 t payload (ballistic coefficient 475), bank 81° with one reversal to the zone's centre-line middle; 5.2 g.
- **Final approach overlay** (`approach_overlay` in `flt15_descents.py`): the last 180 km re-shaped through (26.20 N, 97.20 W) north of the SPI town and (26.10 N, 97.03 W) over the Gulf to the tower, with unit tangents on a chord-length Hermite; the simulated distance is stretched ×1.13. Heights, speeds, times unchanged. Do not use the raw simulated path's last minutes for the page: it reaches the tower straight from the west-north-west.
- **Booster:** ballistic coefficient 3400 kg/m² and the landing-burn start height 4.5 km are assumptions; the landing-burn start point east of the tower is iterated so the booster arrives from the Gulf. Boostback 46 m/s², 1.66 km/s. Cartoon-level fidelity; MECO state from the shaped ascent (66 km, 1.23 km/s ground, 0.9 km/s climb at T+2:19).
- **Zone windows** are the Flight 14 windows re-timed by keeping the offsets from the window open and close (`flt15_zones.retime`); the page and maps say they are not published. Replace `navwarning/` and the windows when Flight 15 warnings appear; the zones are read through `ift14_fit`, so new KMLs need the `NAV` path or a local reader.
- **Template changes made for this mission** (in `~/launch/launch-page-template`, tests pass): `meta.orbit_end_met` (flight ends with a descent), `meta.companion` (second vehicle on its own track, used for Super Heavy), `meta.layer_labels`, `meta.zone_noun`, and segment `kind: 'starship'` in `vehicle-art.js` (the Flight 14 cartoon: nose cone, heat-shield side, flaps, engine skirt) with a sunlit glint in `app.js`. The booster stays a plain stage (user: fine). Documented in `PROFILE_FORMAT.md`. After a template change, copy `docs/js/app.js` and `vehicle-art.js` into this project's `docs/js` or rebuild with `new_page.py`.
- **Maps:** 13.2 in at 300 dpi (3960 px); `ift14_style` fonts scaled 0.6. Sun 0° yellow and Sun −12° cyan, both solid (`linestyles='solid'`, or matplotlib dashes negative levels). Orbit is steel blue so it does not read as the −12° line. The approach map's inset must have `set_zorder(20)` or the parent's city labels show through it.

## Page

`docs/` is a build of the template: do not edit `docs/js` here; change the template, rebuild with `new_page.py`, then `check_page.py`. Quick check passed 9/9 and the full check passed (see work/check_full.log). The Chrome extension was not connected during the build, so the visual check used headless Chrome screenshots (`.shots/`, not published).

## After the launch

Fetch the earliest TLEs with the space-track skill, compare the plane and period with the model (`flt15_tle.txt`), solve T0, replace the zones with the published warnings, and re-run from step 2 with the real deorbit and catch times if they differ.
