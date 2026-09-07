# Coma Virtual Observatory — open items

Live page: https://ulli1991.github.io/coma_flight/ (built from `index.html` + `data/`).
Push to `main` and GitHub Pages redeploys in about 80 s.

## Done on raven (Sept 2026) — regenerate with

```
for s in 139 27 56 77 109 121; do sbatch --export=ALL,SNAP=$s -J cubes_$s run_cubes.sbatch; done
python pack_cubes.py check          # centre + axis order of raw_139.h5 vs the shipped rho384
python pack_cubes.py calib          # once: the check, then fits the z=0 scalings -> cube_scales.json
python pack_cubes.py 139            # dm384 / idm192 / shock384 / ishock192
for s in 27 56 77 109 121; do python pack_cubes.py $s; done   # ep{S}_*192.u8.gz
```

`extract_cubes.py` needs `/ptmp/uli/coma_cubes/{cen139,tracers139}.npy` (BCG centre =
FOF group 0, and the 3000 most bound DM particles of the z=0 BCG that locate the main
progenitor at earlier times). One snapshot takes about a minute on a raven node.
`pack_cubes.py` (and `check`) need a few GB: run them on an interactive node, e.g.
`srun -p interactive -n1 -c8 --mem=24G -t 10 python pack_cubes.py check`.
`extract_dm.py` is gone (it needed freeze_139.h5 from viper, which raven cannot see).

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
- Time evolution (keys , and . or the z strip top right): z = 1.95, 1.0, 0.5, 0.2, 0.1
  from snapshots 27, 56, 77, 109, 121; 192^3 cubes of gas density, star light + colour,
  X-ray, temperature, DM and shocks around the main progenitor, on the z=0 intensity
  scales. Loaded on first use (~20 MB per epoch). Star sprites, cold-gas sprites, black
  holes and labels are z=0 objects and hide at other epochs; the camera stays put.

## Page-side ideas not done yet
- Galaxy labels: NGC 4874 and the NGC 4839 group are placed at the 2nd and 4th
  most massive black holes (400 kpc/h and 2.4 Mpc/h from the core). Replace with
  the real counterparts if known.
- Tour pacing, label placement, ICL stretch, scale-bar width, shock colours and
  brightness, epoch intensity at z = 1.95, the z = 0.1 -> 0 jump: all untested by
  eye (WebGL2 does not work in headless Firefox on raven, so the page can only be
  checked in a real browser).
- The outskirts carry Mach 10-300 accretion shocks that saturate the Mach channel;
  if the outer shell is too loud, raise the upper Mach scale in `cube_scales.json`
  (`mach`, 1..5 now) and repack, or add a radial fade in the shader's mode-5 branch.
- Epochs use the 192^3 grids only (no inner high-res core, no per-epoch sprites);
  a 384^3 density for z = 0.2 / 0.1 would cost ~20 MB each.
- Shock mode shows every shock-finder cell, including Mach 100+ shocks inside galaxy
  ISM (saturated white-yellow dots); a temperature mask would keep only ICM shocks.
- More epochs (z = 3, 1.5, 0.7, 0.35) are one sbatch + pack run each.
- `_bench_baseline.html` / `_bench_expC.html`: untracked benchmark copies of an
  earlier session, safe to delete.

## Controls (current)
drag / arrows look · WASD fly · SPACE/C rise/sink · SHIFT boost · Q/E roll
1 density · 2 X-ray · 3 temperature · 0 dark matter · 4 shocks · I intracluster light
, . step through epochs · 5–9 sights · T tour · H hide HUD · G volume 1:1 · B bloom on/off
P save PNG · L copy link (includes the epoch) · R reset

Photo recipe: frame the view, release all keys, wait for "still N" in the HUD to
pass ~30, press P. Press L to copy a link that reopens the same view.
