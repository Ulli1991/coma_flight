# Coma Flight

**▶ Open the Coma Virtual Observatory: https://ulli1991.github.io/coma_flight/**

The full-resolution build (`index.html`, streams ~30 MB of data from `data/`)
runs on GitHub Pages at the link above: guided tour, gas / X-ray / temperature /
dark-matter / shock / entropy / radio-relic / intracluster-light modes with colour
bars and a probe readout, nine epochs from z = 1.95 to today (keys , and .) with the
galaxies, cold gas, black holes and tracked named galaxies of each epoch, shareable
views (L), PNG export (P), video capture (V).
See `DEEP.md` for what is in it and `TODO.md` for how the cubes are made
(`extract_cubes.py`, `run_cubes.sbatch`, `pack_cubes.py` on raven). The rest of
this file describes the light self-contained artifact version.

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
  T autopilot, H photo mode, adaptive render quality

## Controls
drag / arrows look · WASD fly · SPACE/C rise/sink · SHIFT boost ·
Q/E roll · 1/2/3 density / X-ray / temperature · 5–9 sights · T autopilot ·
H hide HUD · , . step epoch · Z time-lapse z = 2 → 0 · X hide/show galaxies · N hide/show dense gas · P photo (PNG) · V record video · L copy link to this view · I intracluster light ·
0 dark matter · 4 shocks · K entropy · F radio relics · , . step through the epochs (these only in the
GitHub Pages build, and only when their cubes are present) · B bloom on/off · R reset

## Data provenance
Built from /ptmp/uli/movie011/grids/freeze_139.h5, epoch_139{,_vel}.h5,
game_data.npz, game_gas.npz, inner_rho.u8.gz on viper.
`extract_gas.py` and `grid_inner.py` are the extraction jobs (need the
raven mount, run on interactive nodes).

The dark-matter, shock (Arepo shock-finder Mach number + energy dissipation)
and earlier-epoch cubes come from the step_011 snapshots on raven:
`extract_cubes.py` (one Slurm job per snapshot, `run_cubes.sbatch`, about two
minutes each) writes raw float cubes to `/ptmp/uli/coma_cubes/raw_SNAP.h5`, and
`pack_cubes.py` turns them into the `data/*.u8.gz` files on the intensity
scales of the shipped z=0 cubes (`cube_scales.json`). The centre is FoF group 0
of snapshot 139; `pack_cubes.py check` confirms it against the shipped
`data/rho384.u8.gz` by cross-correlation (zero-voxel shift, axis order x,y,z).
Earlier epochs (snapshots 27, 41, 56, 65, 77, 91, 109, 121 = z 1.95, 1.5, 1.0, 0.7,
0.5, 0.35, 0.2, 0.1) follow the main progenitor via the 3000 most-bound DM particles
of the z=0 BCG; each ships its own galaxy sprites (`ep*_stars.bin.gz`).
