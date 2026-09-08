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
         'idm192.u8.gz': 192 ** 3, 'ishock192.u8.gz': 2 * 192 ** 3, 'mag192.u8.gz': 192 ** 3, 'stars.bin.gz': None, 'gas.bin.gz': None}
for s, hi in snaps:
    if s == 139: continue
    p = 'ep%03d_' % s
    files[p + ('pk384.u8.gz' if hi else 'pk192.u8.gz')] = 3 * (384 if hi else 192) ** 3
    for k in ('xray', 'temp', 'dm', 'mag'): files[p + k + '192.u8.gz'] = 192 ** 3
    files[p + 'shock192.u8.gz'] = 2 * 192 ** 3; files[p + 'stars.bin.gz'] = None; files[p + 'gas.bin.gz'] = None
    m = json.load(open(D + '/data/' + p + 'meta.json')); print('%-22s black holes %d, labels %s' % (p + 'meta.json', len(m.get('bh', [])) // 4, [(l['n'], l['p'], l['n_tracers']) for l in m.get('labels', [])]))
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
print('epochs in page:', snaps); print('total data/: %.0f MB, referenced %.0f MB, problems: %d' % (sum(os.path.getsize(D + '/data/' + f) for f in os.listdir(D + '/data')) / 1e6, tot / 1e6, bad))
