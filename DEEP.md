# Coma Virtual Observatory

High-resolution build of the Coma_011 z=0 fly-through, designed for
GitHub Pages (no artifact size cap): `index.html` fetches its data files
from `data/` at runtime — gas at native 384^3 (+ 192^3 inner core),
full-resolution star-light cubes, millions of particles, CGM layer.

Light 16 MB artifact version lives in ~/coma_flight.

Added Sept 2026 from the snapshots on raven (`extract_cubes.py` -> `pack_cubes.py`,
scalings in `cube_scales.json`, all cubes C-ordered x,y,z uint8 over +-4500 ckpc/h):
- `dm384.u8.gz`, `idm192.u8.gz`: dark-matter density (mode 3, key 0)
- `shock384.u8.gz`, `ishock192.u8.gz`: RG cubes, R = Arepo shock-finder Mach number
  (1..5, max in voxel), G = log shock energy dissipation (mode 5, key 4); the shock
  finder carries no shocks in star-forming cells, so the shader shows the accretion
  shock (Mach 10-300 at the edge) and the ICM merger shocks without a temperature mask
- `ep{027,041,056,065,077,091,109,121}_{pk,xray,temp,dm,shock}192.u8.gz`: z = 1.95, 1.5,
  1.0, 0.7, 0.5, 0.35, 0.2, 0.1 around the main progenitor (3000 most bound DM particles
  of the z=0 BCG, snapped to the nearest massive FOF group); pk = density + star light +
  g-r as at z=0; `ep{109,121}_pk384.u8.gz` the same at 384^3 for the two late epochs
- `ep{S}_stars.bin.gz`: the galaxy sprites of each epoch in the `stars.bin.gz` record
  format, encoded to match the shipped z=0 file (`pack_cubes.py calib_sp`, see TODO.md)
