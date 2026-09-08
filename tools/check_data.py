# Integrity check of every file index.html references: decompressed voxel counts, sprite record sanity,
# meta.json contents.  Run after packing (a few GB of RAM):
#   srun -p interactive -n1 -c4 --mem=24G -t 15 /u/uli/jupyter_env/bin/python3 tools/check_data.py
# Integrity check of every file the page references: decompressed sizes, sprite record sanity.
import re, os, gzip, json, numpy as np
D = '/raven/u/uli/coma_flight'
html = open(D + '/index.html').read()
ep = re.search(r'const EPOCHS=\[(.*?)\];', html).group(1)
snaps = [(int(s), 'hi:true' in body) for s, body in re.findall(r'\{s:(\d+),([^}]*)\}', ep)]
files = {'rho384.u8.gz': 384 ** 3, 'slum384.u8.gz': 384 ** 3, 'scol384.u8.gz': 384 ** 3, 'dm384.u8.gz': 384 ** 3, 'shock384.u8.gz': 2 * 384 ** 3,
         'xray192.u8.gz': 192 ** 3, 'temp192.u8.gz': 192 ** 3, 'irho192.u8.gz': 192 ** 3, 'islum192.u8.gz': 192 ** 3, 'iscol192.u8.gz': 192 ** 3,
         'idm192.u8.gz': 192 ** 3, 'ishock192.u8.gz': 2 * 192 ** 3, 'mag192.u8.gz': 192 ** 3, 'stars.bin.gz': None, 'gas.bin.gz': None,
         'lite_pk192.u8.gz': 3 * 192 ** 3, 'lite_dm192.u8.gz': 192 ** 3, 'lite_shock192.u8.gz': 2 * 192 ** 3,
         'vel192.u8.gz': 4 * 192 ** 3, 'met192.u8.gz': 192 ** 3, 'bvec192.u8.gz': 3 * 192 ** 3, 'outer192.u8.gz': 2 * 192 ** 3}
for s, hi in snaps:
    if s == 139: continue
    p = 'ep%03d_' % s
    files[p + ('pk384.u8.gz' if hi else 'pk192.u8.gz')] = 3 * (384 if hi else 192) ** 3
    for k in ('xray', 'temp', 'dm', 'mag', 'met'): files[p + k + '192.u8.gz'] = 192 ** 3
    files[p + 'vel96.u8.gz'] = 4 * 96 ** 3; files[p + 'bvec192.u8.gz'] = 3 * 192 ** 3
    files[p + 'shock192.u8.gz'] = 2 * 192 ** 3; files[p + 'stars.bin.gz'] = None; files[p + 'gas.bin.gz'] = None
    m = json.load(open(D + '/data/' + p + 'meta.json')); print('%-22s black holes %d, labels %s' % (p + 'meta.json', len(m.get('bh', [])) // 4, [(l['n'], l['p'], l['n_tracers']) for l in m.get('labels', [])]))
# the time slider: every snapshot listed in movie.json as a 128^3 RGB cube; the galaxy table
if os.path.exists(D + '/data/movie.json'):
    mv = json.load(open(D + '/data/movie.json')); print('movie.json: %d snapshots, z %.2f .. %.2f, labels in %d' % (len(mv), mv[0]['z'], mv[-1]['z'], sum('lab' in e for e in mv)))
    for e in mv:
        files['mv/mv_%03d.u8.gz' % e['s']] = 3 * 128 ** 3
        if e.get('flow'): files['mv/dv_%03d.u8.gz' % e['s']] = 3 * 48 ** 3
    if os.path.exists(D + '/data/movie_gal.json'):
        mg = json.load(open(D + '/data/movie_gal.json')); print('movie_gal.json: %d snapshots, galaxies per snapshot %d .. %d' % (len(mg), min(len(g) for g in mg), max(len(g) for g in mg)))
if os.path.exists(D + '/data/galaxies.json'):
    gal = json.load(open(D + '/data/galaxies.json')); print('galaxies.json: %d subhaloes, brightest %s' % (len(gal), gal[0]))
if os.path.exists(D + '/data/galaxy_tracks.json'):
    tk = json.load(open(D + '/data/galaxy_tracks.json')); print('galaxy_tracks.json: %d galaxies, median track %d snapshots' % (len(tk), sorted(len(q['t']) for q in tk)[len(tk) // 2]))
if os.path.exists(D + '/data/sky.json'):
    sk = json.load(open(D + '/data/sky.json')); print('sky.json: observer frame of step %d, D = %.0f ckpc/h, %d 2M++ galaxies' % (sk['step'], sk['D'], len(sk['gal'])))
bad = 0; tot = 0
for fn, n in sorted(files.items()):
    path = D + '/data/' + fn
    if not os.path.exists(path): print('MISSING', fn); bad += 1; continue
    raw = gzip.decompress(open(path, 'rb').read()); tot += os.path.getsize(path)
    if n is not None:
        ok = len(raw) == n; bad += not ok
        print('%-22s %6.1f MB gz  %s' % (fn, os.path.getsize(path) / 1e6, 'ok' if ok else 'BAD size %d != %d' % (len(raw), n)))
    else:
        st = fn.endswith('stars.bin.gz'); rec = 10 if st else 8
        ok = len(raw) % rec == 0; bad += not ok; a = np.frombuffer(raw, np.uint8).reshape(-1, rec)
        xyz = a[:, :6].copy().view('<i2').reshape(len(a), 3).astype(np.float64) / 32767.
        r = np.linalg.norm(xyz, axis=1)
        # sprite light within 0.25 (1125 ckpc/h) of the centre: the progenitor should dominate
        w = 0.08 + 2.2 * (a[:, 6] / 255.) ** 2.5 if st else np.ones(len(a))
        print('%-22s %6.1f MB gz  %s  n=%d  median|pos|=%.3f  light r<0.25: %.2f  lum pct 10/50/90: %s  g-r pct 10/50/90: %s' % (
            fn, os.path.getsize(path) / 1e6, 'ok' if ok else 'BAD record', len(a), np.median(r), (w * (r < .25)).sum() / w.sum(),
            np.percentile(a[:, 6], [10, 50, 90]), np.percentile(a[:, 7], [10, 50, 90]) if st else '-'))
allsz = sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(D + '/data') for f in fs)
print('epochs in page:', snaps); print('total data/: %.0f MB, referenced %.0f MB, problems: %d' % (allsz / 1e6, tot / 1e6, bad))
