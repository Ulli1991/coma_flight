# Coma Virtual Observatory — open items

Live page: https://ulli1991.github.io/coma_flight/ (built from `index.html` + `data/`).
Push to `main` and GitHub Pages redeploys in about 80 s.

## Needs the snapshot on raven

### Dark matter (mode 3, key 0)
- Run `extract_dm.py` on raven (needs `/raven/ptmp/uli/sims/.../snapdir_139` and
  `/ptmp/uli/movie011/grids/freeze_139.h5` for the centre and half-width).
- It writes `/ptmp/uli/movie011/dm384.u8.gz` and `idm192.u8.gz`.
- Copy them to `data/dm384.u8.gz` and `data/idm192.u8.gz`, commit, push.
  The page detects them and shows the "0 DM" button; nothing else to change.

### Shocks (planned mode 5, key S)
Two options, in order of quality:
1. Mach-number cube: check whether `PartType0/Machnumber` exists in `snap_139.*.hdf5`
   (Arepo shock finder). If yes, grid it like `extract_dm.py` (384^3 + 192^3 inner,
   max-in-cell rather than mass sum), scale 1..~5 to 0..255, ship as
   `data/mach384.u8.gz` / `data/imach192.u8.gz`.
2. No Mach field: shader-side |grad log(rho*T)| from the cubes already loaded
   (temperature is 192^3, ~47 kpc/h voxels, so only Mpc-scale shocks are crisp),
   or velocity divergence from `epoch_139_vel.h5` on viper.

### Time evolution (big one)
- Extract density + temperature (192^3 is enough, ~7 MB each gzipped) for a handful
  of earlier snapshots of the step_011 film, e.g. z = 2, 1, 0.5, 0.2, 0.
- Page side: a key or slider to swap the cubes, camera stays put.

## Page-side ideas not done yet
- Galaxy labels: NGC 4874 and the NGC 4839 group are placed at the 2nd and 4th
  most massive black holes (400 kpc/h and 2.4 Mpc/h from the core). Replace with
  the real counterparts if known.
- Tour pacing, label placement, ICL stretch, scale-bar width: all untested by eye.
- B key: currently bypasses the whole post stage (no tone map, twice as bright).
  Either remove it or make it a pure bloom toggle.
- `_bench_baseline.html` / `_bench_expC.html`: untracked benchmark copies of an
  earlier session, safe to delete.

## Controls (current)
drag / arrows look · WASD fly · SPACE/C rise/sink · SHIFT boost · Q/E roll
1 density · 2 X-ray · 3 temperature · 0 dark matter (if cubes present) · I intracluster light
5–9 sights · T tour · H hide HUD · G volume 1:1 · B post bypass · P save PNG · L copy link · R reset

Photo recipe: frame the view, release all keys, wait for "still N" in the HUD to
pass ~30, press P. Press L to copy a link that reopens the same view.
