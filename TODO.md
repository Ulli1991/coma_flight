# Coma Virtual Observatory — open items

Live page: https://ulli1991.github.io/coma_flight/ (built from `index.html` + `data/`).
Push to `main` and GitHub Pages redeploys in about 80 s.

## Done on raven (Sept 2026) — regenerate with

```
for s in 139 27 41 56 65 77 91 109 121; do sbatch --export=ALL,SNAP=$s -J cubes_$s run_cubes.sbatch; done
python pack_cubes.py check          # centre + axis order of raw_139.h5 vs the shipped rho384
python pack_cubes.py calib          # once: the check, then fits the z=0 scalings -> cube_scales.json
python pack_cubes.py calib_sp       # once: the galaxy-sprite encoding of stars.bin.gz -> cube_scales.json
python pack_cubes.py calib_cg       # once: the cold-gas sprite encoding of gas.bin.gz -> cube_scales.json
                                    # (label tracking needs /ptmp/uli/coma_cubes/label_tracers.npz, made by
                                    #  tools/label_tracers.py from the z=0 group catalogue; parts run: gas dm stars sp cold bh lab)
python pack_cubes.py 139            # dm384 / idm192 / shock384 / ishock192
for s in 27 41 56 65 77 91 109 121; do python pack_cubes.py $s; done   # ep{S}_*192.u8.gz, ep{S}_stars.bin.gz, ep{109,121}_pk384
```

`extract_cubes.py` needs `/ptmp/uli/coma_cubes/{cen139,tracers139}.npy` (BCG centre =
FOF group 0, and the 3000 most bound DM particles of the z=0 BCG that locate the main
progenitor at earlier times). One snapshot takes about a minute on a raven node.
Parts can be added to an existing raw file without redoing the rest:
`sbatch --export=ALL,SNAP=109,PARTS=sp+hires -J cubes_109 run_cubes.sbatch`
(parts: gas dm stars sp hires; join them with `+`, Slurm eats commas).
`pack_cubes.py` (and `check`) need a few GB: run them on an interactive node, e.g.
`srun -p interactive -n1 -c8 --mem=24G -t 10 python pack_cubes.py check`.
`extract_dm.py` is gone (it needed freeze_139.h5 from viper, which raven cannot see).

Sept 8: galaxy sprites per epoch, three more epochs, 384^3 density for z = 0.2 / 0.1, colour bars.
- `ep{S}_stars.bin.gz`: the star particles of that snapshot in the `stars.bin.gz` record
  (int16 xyz / 32767 * 4500 ckpc/h from the epoch's centre, bytes lum, g-r, h=90, pad).
  The shipped z=0 file was made on viper from data raven does not have; its lum / g-r
  bytes are noisy, non-linear functions of the particle's own r-band L and g-r (r ~ 0.7
  per particle, no neighbourhood sum reproduces them either). `pack_cubes.py calib_sp`
  therefore matches distributions: threshold = the 939637-th brightest particle at z=0
  (log L_r = 4.00, mag zero point), and quantile-quantile lookup tables byte(log L_r),
  byte(g-r) (`sp_lum`, `sp_col` in `cube_scales.json`) so the regenerated z=0 file has
  exactly the shipped byte histograms (total sprite light 0.999 of shipped; 96^3 light
  maps correlate at 0.68 because the shipped selection is noisy). Young populations are
  much brighter, so every epoch before z = 0.1 has more particles above the z=0 threshold
  than the cap of 1.2M (`SP_CAP`): those files hold the 1.2M brightest (cut at log L_r
  4.64 at z = 1.95 ... 4.18 at z = 0.35), 7-8 MB each. `/ptmp/uli/coma_cubes/stars139_regen.bin.gz`
  is the regenerated z=0 file (not shipped): copy it over `data/stars.bin.gz` to have all
  epochs on the identical pipeline.
- Epochs z = 1.5, 0.7, 0.35 (snapshots 41, 65, 91). z = 3 does not exist: the run's
  first snapshot is 27 at z = 1.95. Lookback times from h = 0.681, Om = 0.306.
- `ep109_pk384.u8.gz`, `ep121_pk384.u8.gz`: density + star light + g-r at 384^3
  (`hi: true` in `EPOCHS`); X-ray, temperature, DM and shocks stay 192^3 as at z=0.
- Shocks made visible (the request: the cluster's accretion shock and the ICM merger shocks).
  Checked on the z=0 cells: the Arepo shock finder reports no shock in star-forming
  cells (extract_cubes.py now guards against it anyway, a no-op here), and of the 42k
  cells with Mach >= 10 all but 22 lie beyond a third of the cube at 1e5-1e7 K: that
  is the accretion shock. The shader's old temperature mask (nothing below ~2e7 K)
  therefore hid the whole accretion shock and most shocks before z = 0.5; the ICM
  merger shocks (Mach 2-5 inside, 1e7-6e7 K) survive without it. Shader: colour ramp Mach 1.2 .. 4
  (`smoothstep(0.05,0.75,s.r)`), emissivity `x*(0.4+1.6*g)*6.5` instead of
  `x*(0.25+1.75*g^2)*2.2` (a sheet is 1-2 ray steps thick where the gas fills 100), no
  radial dimming of the accretion shell any more, and an exposure floor of 0.4 in shock
  mode (`expoV`) because the exposure falls to ~0.1 near the core where the merger
  shocks are. `cube_scales.json` `ed` was refitted after the ISM removal. All untested by
  eye: if the accretion shell now swamps the core, bring back
  `x*=1.-0.6*smoothstep(0.5,0.9,redge)`; if the sheets are too loud, lower 6.5.
- Colour bar (bottom right, hidden with the HUD): the shader ramp of the current field with
  ticks in physical units through the same smoothstep stretches: n_H [cm^-3] for the gas
  (voxel mass -> density with h = 0.681, X = 0.76; comoving voxels, so physical only at
  z=0), relative X-ray emissivity per dex, T [K], DM density [Msun/kpc^3], g-r for the
  ICL, Mach number for the shocks. Constants copied from `cube_scales.json` into `SCL`.

Two parallel sessions on 2026-09-07 built this twice; the other attempt
(`extract_raven.py`, branch `worktree-raven-todo`, outputs in
`/raven/ptmp/uli/coma_flight_cubes/`) is superseded by this one: it had only density
and temperature at the earlier epochs and a single-channel Mach cube. Its centre
cross-check survives as `pack_cubes.py check`; its page used keys `[` `]` and the hash
key `e` for the epoch, this page uses `,` `.` and `z`.

- Dark matter (mode 3, key 0): `data/dm384.u8.gz`, `data/idm192.u8.gz`.
- Shocks (mode 5, key 4): Arepo shock finder `PartType0/Machnumber` (max per voxel) and
  `EnergyDissipation` (sum) packed as an RG cube, `data/shock384.u8.gz` + `data/ishock192.u8.gz`.
  Colour = Mach number (1..5), brightness = dissipation, faint grey gas for orientation.
- Time evolution (keys , and . or the z strip top right): z = 1.95, 1.5, 1.0, 0.7, 0.5,
  0.35, 0.2, 0.1 from snapshots 27, 41, 56, 65, 77, 91, 109, 121; 192^3 cubes of gas
  density, star light + colour, X-ray, temperature, DM and shocks around the main
  progenitor (384^3 density + star light at z = 0.2, 0.1), on the z=0 intensity scales,
  plus that snapshot's galaxy sprites, cold-gas sprites, black holes and the tracked
  named galaxies (ep{S}_meta.json). Loaded on first use (~30 MB per epoch, ~50 MB for
  the two 384 ones); the sights (J / label clicks) are z=0 only; the camera stays put.

Sept 8, second round (page: entropy, radio, probe, legend, video, texture freeing; raven: cold gas,
black holes and tracked labels per epoch):
- Entropy mode (key K, mode 6): K = kT / n_e^(2/3) [keV cm^2] computed in the shader from the
  temperature and density bytes (`uKt` constants from `cube_scales.json`), shown over
  200..6300 keV cm^2, dark where the gas is below the temperature scale. No data needed.
- Radio relic mode (key F, mode 7, needs the shock cubes): emissivity = dissipation (linear
  again, the byte is log over ~3 dex) x efficiency(Mach) with the efficiency a smoothstep
  Mach 1.8..3.6 (Hoeft & Bruggen 2007 in spirit). Merger shocks in dense gas light up, the
  accretion shock does not. Both modes have colour bars and tour legs; `#m=6/7` in links.
- HUD probe ("here"): n_H, T and rho_DM at the camera position from the outer cubes kept
  in memory (`PROBE0`, `q.arr`; ~130 MB for z=0, ~35 MB per resident epoch).
- Legend under the colour bar for the galaxy colours, cold gas and black holes.
- V records the canvas to WebM (MediaRecorder, 30 fps, 24 Mbit/s); V again saves the file.
- Only the current epoch and its neighbours stay on the GPU (`freeEpoch`); a freed epoch is
  re-fetched from the browser cache when stepped to again.
- `ep{S}_gas.bin.gz`: the 200k densest cold / star-forming cells of the snapshot in the
  `gas.bin.gz` record (int16 xyz, byte log rho via a quantile table `cg_lrho` fitted to the
  shipped file, byte SF flag 0/255). The shipped byte is linear in log rho (r = 0.98).
- `ep{S}_meta.json`: `bh` = [x, y, z, m, ...] of every black hole in the cube, m as in `D_BH`
  (m = (log10 M_BH[Msun/h] - 5.833) / 4.419, floor 0.01583 = seed mass; fit r = 1.000), and
  `labels` = the tracked positions of NGC 4889, NGC 4874 and the NGC 4839 group: the
  median of the 1000 most bound DM particles of their z=0 subhaloes (`label_tracers.npz`,
  subhaloes 0, 75, 1 of group 0) found in that snapshot (`extract_cubes.py SNAP lab`).
  The page moves the labels with them; the z=0 label positions are now these subhalo
  centres (NGC 4839 moved 136 kpc/h from the black-hole position used before).
- Untested by eye: entropy and radio stretches, the probe line's length in the HUD, the
  legend row, whether the tracked labels sit on the right galaxies at early epochs
  (NGC 4874 is a 7e10 Msun stripped satellite today; its tracers may scatter early on,
  `n_tracers` in meta.json tells how many were found), video capture on Safari (mp4).

Sept 8, third round (magnetic field, Compton y, one key scheme, labels, mode crossfade):
- Magnetic field (key 7, mode 8): `data/mag192.u8.gz` and `ep{S}_mag192.u8.gz` = log10 of the
  mass-weighted |B| [uG] of the non-star-forming gas per 192^3 voxel (`extract_cubes.py SNAP mag`,
  `pack_cubes.py mag`; Arepo code units -> Gauss with 2.60e-6 h / a^2, h = 0.681), physical uG at
  every epoch over 0.001 .. 10 uG (`mag` in `cube_scales.json`; the floor was 0.01 uG until the z = 2 audit showed only 6% of the voxels above it). z=0 profile: 3.3 uG inside 100 kpc/h,
  1.2 at 300-600, 0.34 at 1-1.5 Mpc/h, 0.11 at 1.5-2.5, 0.017 beyond, peak 11 uG. Shader: colour by
  field strength (`bmap`, indigo -> violet -> magenta -> pale yellow), emissivity = gas density x
  (0.3 + 1.7 x). The star-forming ISM is left out because its 10s of uG would mark every galaxy.
  The HUD probe shows B. Loaded like DM / shocks: the button appears when the cube loads; an epoch
  without its cube falls back to density.
- Compton y (key 5, mode 9): no data, computed in the shader from the density and temperature bytes as
  the electron pressure P_e = 1.17 n_H kT (`PY` constants, `uPy`), shown over 3 dex below the top of the
  byte ranges (~3e-6 .. 4e-3 keV cm^-3 on the 23 kpc/h voxel means); dark below the temperature scale
  so the cold ISM does not show its density. Colour `ymap` (black -> blue -> sky -> white).
- One key scheme: the fields sit on the digit keys in the order 1 DM, 2 density, 3 temperature,
  4 X-ray, 5 Compton y, 6 entropy, 7 magnetic, 8 shocks, 9 radio, 0 ICL (`ORDER`; the internal mode ids
  and the `#m=` of old links are unchanged). The letter aliases I / K / F are gone. The sights moved
  from 5-9 to J (cycles through them) and to clicking a label.
- Labels no longer overlap: named galaxies first (nearest first), then the sights; a label that would
  cover a placed label or a HUD panel climbs up its leader line (up to 150 px, then hides).
- Mode switches crossfade (0.5 s) through the same hold buffer as the epoch switches (0.9 s).
- Untested by eye: the magnetic stretch (if the core is a flat magenta blob lower the 1.7 boost or the
  `hi`), the Compton-y range (`PY[3]`, 3 dex), the label climbing, the ten buttons on a narrow window.

Sept 8, fourth round ("beyond state of the art"): new physics, the cosmic web, the time slider, the inspector.
- Kinematics (`extract_cubes.py SNAP kin`, `pack_cubes.py kin` -> `vel192.u8.gz` at z = 0, `ep*_vel96.u8.gz` at the epochs (the field varies over the whole box, 19 MB at 192^3), RGBA):
  mass-weighted velocity of the hot (T > 1e6 K, non-SF) gas relative to the systemic velocity (the hot
  gas within 500 kpc/h), byte 128 + v/1500 km/s x 127, and the in-voxel dispersion sigma = sqrt(<v^2> - <v>^2)
  as log10 over 10 .. 1000 km/s. Modes: K velocity along the line of sight (blue approaching, +-800 km/s,
  colour bar), U turbulence (sigma). The XRISM Coma result (~200 km/s bulk, ~ 200 km/s dispersion) is the
  comparison. The probe shows sigma.
- Metals (`met`): mass-weighted Z/Zsun (Zsun = 0.0127) of the non-SF gas, log10 over -2 .. 0.5, mode M.
- B vector (`bvec`, RGB, signed asinh encoding 128 + 127 asinh(B/0.05 uG)/asinh(200)): F = Faraday rotation
  measure, a screen integral in the shader, RM = 812 sum n_e B_par dl [rad/m^2] with n_e from the density
  byte where hot gas exists and dl = 4500/0.681/(1+z) kpc per unit; coloured after the ray loop, blue
  negative / red positive over +-300 (RMLIM), white beyond. O = radio halo, emissivity ~ n sigma^2
  B^2/(B^2 + B_CMB^2) with B_CMB = 3.24 (1+z)^2 uG (turbulent re-acceleration in spirit).
- The cosmic web (`extract_cubes.py 139 outer`, `pack_cubes.py outer` -> `outer192.u8.gz`, RG = log total
  matter, log gas over +-30 Mpc/h, 312 kpc/h voxels, all particle types incl. the low-res ones): the
  shader's wide branch takes over when the camera is outside r = 1.02 (the flight speed scales with r);
  inside r < 0.97 the fine density cube is drawn. Tour leg 2. HUD says "cosmic web".
- Time slider (MOVIE chip in the epoch strip, or Z): `extract_movie.py` (job array `run_movie.sbatch`,
  or the interactive partition when the group CPU limit blocks the array) makes 128^3 cubes of density,
  log T and star light for every snapshot 27 .. 139 with the same centring and label tracking;
  `pack_cubes.py movie` -> `data/mv/mv_SNAP.u8.gz` (RGB, T and light at 64 levels) + `data/movie.json`
  (z, lookback time, M200, R200, centre offset, label positions). The page streams them around the slider
  (8 ahead, 2 behind, at most 26 resident on the GPU) and blends the two neighbours in the shader;
  only density / temperature / entropy / Compton y exist there (`MVOK`), other modes are refused with a
  caption; the tracked galaxies leave their orbits on an overlay canvas (past solid, future faint).
  Z plays at 3.1 snapshots/s; the tour's assembly leg uses it. The nine full-resolution epochs are untouched.
  Size: 173 MB for the 113 cubes (0.6 MB at z = 2 to 2.3 MB at z = 0); data/ is now 570 MB in all. To halve
  it, pack every other snapshot (`glob` in `pack_movie`, the page needs no change: it reads movie.json).
- Galaxy inspector: click a galaxy at z = 0 (`data/galaxies.json`, `pack_cubes.py galaxies`: the 803
  subhaloes with M* >= 1e9 Msun inside the cube) for its catalogue entry; FLY THERE; Escape closes.
- Radial profiles (I): n_H, T, K, B, sigma, Z, rho_DM spherical means from the resident cubes (computed on
  first use, ~1 s at 384^3), camera radius and R200 marked; works at the epochs too.
- The view from Earth (E, or the EARTH VIEW button): the observer's frame of this realisation comes from
  `~/coma_300/data/zoom_rotations.npz` (step 11: the observer at the parent box centre, MUSIC only
  translates the zoom, North = +x of the sim axes) and `parent_radec.npz` (D = 72154 ckpc/h, the halo
  at RA 193.48 Dec 29.44 vs Coma's 194.95 / 27.98), packed by `pack_cubes.py sky` into `data/sky.json`
  together with the real 2M++ galaxies (Lavaux & Hudson 2011, 6 deg cone, 4000 < v_cmb < 10000 km/s) as
  tangent-plane offsets from NGC 4889. The camera sits at -D n, looks along n with North up and East
  left (screen right = n x up), field of view 0.075 (W / S zoom), the fake sky off, sprites boosted by
  D^2; every field works there (the line-of-sight velocity is then the real Earth line of sight); the
  amber circles are the 2M++ galaxies (size by Ks), with a compass and a 1-degree bar. Known: this
  realisation's massive companion is 1.9 deg east of the core (PA 100), not at NGC 4839's PA 228 /
  0.61 R200; NGC 4874's analogue sits 0.21 deg west. The eROSITA and Planck images of coma_300 are not on
  raven (only the Planck y profile), so no image overlay yet.
- Temperature scale extended (Sept 8, late): the shipped z=0 temp192 started at 1.3e7 K, so at z = 2 (progenitor
  gas ~1e7 K) 99% of the voxels were byte 0 = "cold": the young cluster was cyan-tinted and veiled, and the
  temperature / entropy / Compton-y modes were dark. `pack_cubes.py temp` now writes temp192 for z = 0 (from
  raw_139, replacing the viper cube) and every epoch on log10 T = 5.5 .. 8.0318 (3e5 .. 1.1e8 K; `temp` in
  cube_scales.json, the old scale kept as `temp_viper`; do not re-run `calib`, it would refit the old one), and
  `movie` was repacked. Shader thresholds are now physical: cold tint below ~1e6 K, hot-wind tint 5e7..7e7 K,
  the temperature ramp 3e6 .. 9e7 K, entropy and Compton y only above ~2e6 K. The cyan tint at z = 0 therefore
  marks truly cold gas (1e4 K tails and filaments), no longer everything below 1.7e7 K.
- Floors at high redshift (audit of the z = 1.95 cubes against z = 0): X-ray 8% nonzero (core byte 184), DM 38%, metals
  11%, star light 0.5%, sigma 18% -- the cores sit mid-scale, only the outskirts fall below, physical. Shocks: 56% of
  the z = 2 voxels carry a Mach number vs 10% today (accretion everywhere); if the shock mode is a wall there, raise
  the lower edge of the Mach ramp for the epochs. Every byte scale is calibrated on Coma today: check a new field at
  z = 2 before shipping it.
- Untested by eye: all of it; the velocity colour range (VLIM), RMLIM, the wide-view
  stretches (smoothstep 0.15 / 0.2 on the outer bytes), the movie brightness (the 128^3 voxel is 3x the
  384 one: +3 log10 3 on the mass fields), the movie streaming on a slow connection.

## RAVEN — needs the snapshots (run there, then copy into data/ and push)
- Shock cubes: if Mach 1.2-1.5 turbulence still litters the ICM, raise the lower edge of
  the shader's Mach ramp; the cube scale (`mach`, 1..5) only needs repacking if the
  colour ramp should reach beyond Mach 5.
- If the epoch sprites should keep every particle above the z=0 threshold instead of
  the 1.2M brightest, raise `SP_CAP` in pack_cubes.py (z = 1.95 would need ~3M, 18 MB).

## PAGE — laptop only, needs a browser check (no data involved)
- Galaxy labels: NGC 4874 and the NGC 4839 group sit at the 2nd and 4th most massive
  black holes (400 kpc/h and 2.4 Mpc/h from the core). Move them if the real
  counterparts are known.
- Untested by eye: tour pacing incl. the new assembly leg, label placement, ICL
  stretch, scale-bar width, shock colours after the temperature mask and radial fade,
  epoch intensity at z = 1.95, the crossfade length (0.9 s), time-lapse cadence (3.8 s),
  the epoch sprites (do the z = 0.1 galaxies look like the z = 0 ones across the
  crossfade? if not, swap in stars139_regen.bin.gz, see above), the colour bar
  (placement above the key legend, tick labels near the edges, the n_H tick values).
- The tour's assembly leg downloads all eight epochs (~280 MB) on first play; on a slow
  connection consider skipping that leg until the epochs are cached, or thinning the
  time-lapse to every other epoch. With every epoch loaded the 3D textures take ~1.2 GB
  of GPU memory (two 384^3 pk cubes at 170 MB each); if integrated GPUs choke, drop
  `hi: true` or free the textures of epochs far from the current one.
- `_bench_baseline.html` / `_bench_expC.html`: untracked benchmark copies, safe to delete.

## Controls (current)
drag / arrows look · WASD fly · SPACE/C rise/sink · SHIFT boost
, . step epoch · Z time-lapse · X galaxies on/off · N dense gas on/off
fields on the digit keys: 1 dark matter · 2 density · 3 temperature · 4 X-ray · 5 Compton y · 6 entropy ·
7 magnetic field · 8 shocks · 9 radio relics · 0 intracluster light (the top-right buttons in the same order)
K velocity · U turbulence · M metals · F Faraday rotation · O radio halo · I radial profiles · click a galaxy
MOVIE chip / Z: the 113-snapshot time slider (Z plays / pauses, , . return to the epochs)
J next sight, or click a label · T tour · H hide HUD · G volume 1:1 · B bloom on/off
P save PNG · V record video · L copy link (includes the epoch) · R reset

Photo recipe: frame the view, release all keys, wait for "still N" in the HUD to
pass ~30, press P. Press L to copy a link that reopens the same view.
