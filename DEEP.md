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
  (1..5, max in voxel), G = log shock energy dissipation (mode 5, key 4)
- `ep{027,056,077,109,121}_{pk,xray,temp,dm,shock}192.u8.gz`: z = 1.95, 1.0, 0.5, 0.2,
  0.1 around the main progenitor (3000 most bound DM particles of the z=0 BCG,
  snapped to the nearest massive FOF group); pk = density + star light + g-r as at z=0
