# Flight 15 navigational warnings

The six NGA warnings for Starship Flight 15 as issued on 7–9 Nov 2026 and the Mexican NOTAM for the ship's return, cached
12 Nov (text as broadcast; `work/flt15_zones.py` parses both formats into the polygons used by the page and the maps). The Flight 14 KML polygons that stood in for them until 9 Nov were removed; they
are in `~/sda/starship-ITF-14/navwarning` and in this repository's history.

| file | area | window (UTC) | model use |
|---|---|---|---|
| `NAVAREA-IV-1080-26` | A: Gulf of America, Yucatán Channel, Caribbean south of Cuba; B: off Starbase | 2245–0040Z 13 Nov, alternate daily 14–19 Nov | ascent; Super Heavy return |
| `NAVAREA-IV-1081-26` | Laguna Madre, Padre Island and the Gulf off Starbase (ship catch) | 2108–0014Z 14–15 Nov, alternate daily 15–20 Nov | final approach and catch, T+22:28–22:32:34 |
| `B4288-26` | NOTAM, Mexico FIR (MMFR), SFC–UNL, areas C2 and C3: Pacific off Baja California (120.6°W) across northern Mexico to the Rio Grande north-west of Laredo (99.9°W) | 2056–0014Z daily 14–20 Nov | nominal entry, T+22:19:27–22:26:05 at 72–53 km |
| `HYDROPAC-3244-26` | eastern Indian Ocean | 2253–0114Z daily 13–19 Nov | no-insertion contingency, splash T+56:16 |
| `HYDROPAC-3257-26` | North Pacific, Midway Islands (the Flight 14 North Pacific polygon) | 0130–0501Z daily 14–20 Nov | early Hawaii, orbit 2 (window also covers orbit 3) |
| `HYDROPAC-3254-26` | South Pacific, French Polynesia to west of Chile | 0812–1137Z daily 14–20 Nov | mid-mission Chile, orbit 7 (also orbit 8) |
| `HYDROPAC-3256-26` | North Pacific, Micronesia to north-east of Hawaii | 2217–0144Z daily 14–20 Nov | late Hawaii, orbit 16, the pass after the Starbase pass (also orbit 17) |

The dates are the UTC days on which each window opens; the launch window is 13 Nov 22:45–00:00Z, so the windows dated
14–20 Nov belong to launches on 13–19 Nov. Relative to a 22:45Z launch the contingency windows open at T+02:45 (Midway),
T+09:27 (Chile) and T+23:32 (late Hawaii): in each case within a minute of the moment the modeled descent enters the area at
about 80 km, and they close about one orbit plus 20 minutes after the first-pass splashdown. HYDROPAC 3256/26 therefore
places the late-Hawaii option one orbit after the Starbase pass, not before it. The ship-catch area opens 21:08Z, 9.5 min
before the earliest catch, and stays open one orbit beyond the latest one.

NOTAM B4288/26 opens at 20:56Z, the official entry time (T+22:10:46) after a 22:45Z launch, 8.5 min before the ship
reaches the corridor, and closes with the ship-catch area at 00:14Z. Its areas are numbered C2 and C3 and are cut by the
FIR boundary in the west and the border in the east: a C1 and any notice for the Texas part are not cached here. The
corridor showed that the model's single-bank return ran 50 to 70 km too far south over Mexico; the nominal was re-solved
to follow it (see the main README).
