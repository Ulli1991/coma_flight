# How the shipped cold-gas sprites (data/gas.bin.gz) and black holes (D_BH in index.html) are encoded,
# by matching them to the snapshot-139 cells / PartType5 by position. Result (Sept 2026): the log-rho byte is
# linear in log10 rho (r = 0.98, lo -3.21 hi -0.88, floor 9), the SF byte is 0/255, the shipped set is the
# 200k densest cold or star-forming cells; the BH byte is (log10 M_BH[Msun/h] - 5.833) / 4.419 floored at
# 0.01583 (r = 1.000). pack_cubes.py calib_cg / meta_json use these.  Run on an interactive node:
#   srun -p interactive -n1 -c4 --mem=48G -t 25 /u/uli/jupyter_env/bin/python3 tools/fit_gas_bh_sprites.py
import numpy as np, h5py, glob, gzip, re, json, os, sys, time
from scipy.spatial import cKDTree
OUTD = '/ptmp/uli/coma_cubes'
base = '/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/'
files = sorted(glob.glob(base + 'snapdir_139/snap_139.*.hdf5'), key=lambda p: int(p.split('.')[-2]))
gfiles = sorted(glob.glob(base + 'groups_139/fof_subhalo_tab_139.*.hdf5'), key=lambda p: int(p.split('.')[-2]))
cen = np.load(OUTD + '/cen139.npy'); BOX = 681000.; HALF = 4500.
KB = 1.3807e-16; MP = 1.6726e-24; GAMMA = 5. / 3.
def wrap(c): d = c - cen; d -= np.round(d / BOX) * BOX; return d
html = open('/raven/u/uli/coma_flight/index.html').read()
# ---------------------------------------------------------------- cold gas sprites
raw = gzip.decompress(open('/raven/u/uli/coma_flight/data/gas.bin.gz', 'rb').read()); n = len(raw) // 8
rec = np.frombuffer(raw, np.uint8).reshape(n, 8)
gxyz = rec[:, :6].copy().view('<i2').reshape(n, 3).astype(np.float64) / 32767. * HALF
print('gas.bin.gz: %d records, |pos| pct 50/99/100 (half units): %s' % (n, np.round(np.percentile(np.linalg.norm(gxyz, axis=1) / HALF, [50, 99, 100]), 3)))
print('  byte6 (lrho) pct 0/1/10/50/90/99/100:', np.percentile(rec[:, 6], [0, 1, 10, 50, 90, 99, 100]), ' byte7 (sf) unique:', np.unique(rec[:, 7])[:10], 'frac 255: %.3f' % (rec[:, 7] == 255).mean())
t0 = time.time(); P = []; R = []; S = []; T = []
for fn in files:
    with h5py.File(fn) as f:
        g = f['PartType0']; d = wrap(g['Coordinates'][:].astype(np.float64)); m = np.max(np.abs(d), axis=1) < HALF
        rho = g['Density'][:][m]; u = g['InternalEnergy'][:][m]; xe = g['ElectronAbundance'][:][m]; XH = g['GFM_Metals'][:, 0][m]; sfr = g['StarFormationRate'][:][m]
        mu = 4. / (1. + 3. * XH + 4. * XH * xe); tt = (GAMMA - 1.) * u * 1e10 * mu * MP / KB
        sel = (sfr > 0) | (tt < 10 ** 4.9)
        P.append(d[m][sel]); R.append(rho[sel]); S.append(sfr[sel] > 0); T.append(tt[sel])
P = np.concatenate(P); R = np.concatenate(R); S = np.concatenate(S); T = np.concatenate(T)
print('cold/SF cells in cube: %d (SF %d), %.0f s' % (len(P), S.sum(), time.time() - t0))
tree = cKDTree(P); dist, idx = tree.query(gxyz, distance_upper_bound=0.3); ok = np.isfinite(dist)
u_, cnt = np.unique(idx[ok], return_counts=True); uniq = ok.copy(); uniq[ok] &= ~np.isin(idx[ok], u_[cnt > 1])
j = idx[uniq]; lr = np.log10(R[j]); b6 = rec[uniq, 6].astype(float); b7 = rec[uniq, 7]
print('  matched %d records within 0.3 kpc/h (unique %d); dist pct 50/90: %s' % (ok.sum(), uniq.sum(), np.round(np.percentile(dist[ok], [50, 90]), 3)))
a, b = np.polyfit(lr, b6, 1); r = np.corrcoef(lr, b6)[0, 1]; print('  lrho byte = %.3f log10(rho) + %.3f  r=%.4f rms=%.2f -> lo=%.4f hi=%.4f' % (a, b, r, (b6 - a * lr - b).std(), -b / a, (255 - b) / a))
print('  sf byte vs sfr>0: byte 255 & SF %d, byte 255 & not SF %d, byte 0 & SF %d, byte 0 & not SF %d' % (((b7 == 255) & S[j]).sum(), ((b7 == 255) & ~S[j]).sum(), ((b7 == 0) & S[j]).sum(), ((b7 == 0) & ~S[j]).sum()))
print('  selection: min log rho matched %.3f; cells denser than that: %d (records %d); SF fraction among matched %.3f' % (lr.min(), (np.log10(R) >= lr.min()).sum(), n, S[j].mean()))
QK = [0, 0.1, 0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98, 99, 99.5, 99.9, 99.99, 100]
kx = np.percentile(np.sort(lr), QK); ky = np.percentile(np.sort(b6), QK); print('  q-q knots lrho->byte:', ' '.join('%.2f:%.0f' % (x, y) for x, y in zip(kx, ky)))
# ---------------------------------------------------------------- black holes
D_BH = np.array(json.loads(re.search(r'const D_BH = (\[.*?\]);', html).group(1))).reshape(-1, 4)
bxyz = D_BH[:, :3] * HALF; bm = D_BH[:, 3]
BP = []; BM = []
for fn in files:
    with h5py.File(fn) as f:
        if 'PartType5' not in f: continue
        d = wrap(f['PartType5/Coordinates'][:].astype(np.float64)); m = np.max(np.abs(d), axis=1) < HALF
        BP.append(d[m]); BM.append(f['PartType5/BH_Mass'][:][m])
BP = np.concatenate(BP); BM = np.concatenate(BM)
dist, idx = cKDTree(BP).query(bxyz); print('black holes: page %d, snapshot in cube %d, match dist pct 50/90/100: %s' % (len(bxyz), len(BP), np.round(np.percentile(dist, [50, 90, 100]), 3)))
lm = np.log10(BM[idx] * 1e10); a, b = np.polyfit(lm, bm, 1); print('  m = %.4f log10(M_BH/Msun/h... 1e10 units*1e10) + %.4f  r=%.4f  -> lo=%.3f hi=%.3f; min page m %.5f at logM %.3f; snapshot BH logM pct 0/50/100: %s' % (
    a, b, np.corrcoef(lm, bm)[0, 1], -b / a, (1 - b) / a, bm.min(), lm[np.argmin(bm)], np.round(np.percentile(np.log10(BM * 1e10), [0, 50, 100]), 2)))
print('  page m pct 0/10/50/90/100:', np.round(np.percentile(bm, [0, 10, 50, 90, 100]), 4), ' unmatched snapshot BHs heavier than lightest page BH:', (np.log10(BM * 1e10) > lm.min()).sum() - len(bxyz))
