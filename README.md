# Coma Flight

Interactive WebGL2 free flight through the z=0 snapshot of the Coma_011
constrained re-simulation (TNG100 resolution) — the closing frame of the
step_011 formation film.

`coma_flight.html` is fully self-contained (~15 MB): open it in any
WebGL2 browser, or publish it as a claude.ai artifact.
Live copy: https://claude.ai/code/artifact/04ee4360-12d9-4361-bb20-76108562c6c6

## Contents of the page
- 192^3 gas density cube (+ 192^3 high-res inner core, ±1200 kpc regridded
  from particles), 128^3 X-ray and temperature cubes — gzipped, decoded
  with the native DecompressionStream API
- 260k brightest star particles as point masses (g−r as age proxy:
  blue young → red old), 16-bit positions
- 200k densest cold/star-forming gas cells (cyan layer, jellyfish tails)
- HUD in real units (r_BCG, r/R200, speed), sights on keys 5–9,
  T autopilot, H photo mode, M minimap, adaptive render quality

## Controls
drag / arrows look · WASD fly · SPACE/C rise/sink · SHIFT boost ·
Q/E roll · 1/2/3 density / X-ray / temperature · 5–9 sights · T autopilot ·
H hide HUD · M map · R reset

## Data provenance
Built from /ptmp/uli/movie011/grids/freeze_139.h5, epoch_139{,_vel}.h5,
game_data.npz, game_gas.npz, inner_rho.u8.gz on viper.
`extract_gas.py` and `grid_inner.py` are the extraction jobs (need the
raven mount, run on interactive nodes).
