# tools

Analysis and verification scripts behind the raven pipeline (`extract_cubes.py`,
`pack_cubes.py`). Each file starts with what it found and the `srun` line to rerun it.

- `label_tracers.py` builds `label_tracers.npz` (needed by `extract_cubes.py SNAP lab`)
- `check_data.py` verifies every file the page references after packing
- `fit_star_sprites.py`, `fit_star_sprites_neighbours.py` why the galaxy sprites are
  matched by distribution (`pack_cubes.py calib_sp`)
- `fit_gas_bh_sprites.py` the cold-gas and black-hole encodings (`calib_cg`, `meta_json`)
- `shock_cells.py` where the shocks are (why the shader has no temperature mask)
