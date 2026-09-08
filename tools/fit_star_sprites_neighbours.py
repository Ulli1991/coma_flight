# Follow-up to fit_star_sprites.py: are the shipped lum / g-r bytes sums or means over the k nearest
# star particles or within a radius?  Result (Sept 2026): no, best r = 0.76 for 4 neighbours.
#   srun -p interactive -n1 -c8 --mem=60G -t 25 /u/uli/jupyter_env/bin/python3 tools/fit_star_sprites_neighbours.py
import numpy as np, h5py, glob, os, sys, time, gzip
from scipy.spatial import cKDTree
SCR = os.path.dirname(os.path.abspath(__file__))
base = '/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/snapdir_139/'
files = sorted(glob.glob(base + 'snap_139.*.hdf5'), key=lambda p: int(p.split('.')[-2]))
BOX = 681000.0
cen = np.load('/ptmp/uli/coma_cubes/cen139.npy'); HALF = 4500.
raw = gzip.decompress(open('/raven/u/uli/coma_flight/data/stars.bin.gz', 'rb').read()); n = len(raw) // 10
rec = np.frombuffer(raw, np.uint8).reshape(n, 10)
xyz = rec[:, :6].copy().view(np.int16).reshape(n, 3).astype(np.float64) / 32767. * HALF
lum, gr = rec[:, 6].astype(float), rec[:, 7].astype(float)
P = []; PH = []; A = []; M = []
for fn in files:
    with h5py.File(fn) as f:
        if 'PartType4' not in f: continue
        s = f['PartType4']
        c = s['Coordinates'][:].astype(np.float64); d = c - cen; d -= np.round(d / BOX) * BOX
        m = np.max(np.abs(d), axis=1) < HALF + 50
        if not m.any(): continue
        P.append(d[m]); PH.append(s['GFM_StellarPhotometrics'][:][m]); A.append(s['GFM_StellarFormationTime'][:][m]); M.append(s['Masses'][:][m])
P = np.concatenate(P); PH = np.concatenate(PH); A = np.concatenate(A); M = np.concatenate(M)
L_r = 10 ** (-0.4 * PH[:, 5].astype(np.float64)); L_g = 10 ** (-0.4 * PH[:, 4].astype(np.float64))
tree = cKDTree(P)
dist, idx = tree.query(xyz, distance_upper_bound=0.3); ok = np.isfinite(dist)
u, cnt = np.unique(idx[ok], return_counts=True); shared = set(u[cnt > 1])
uniq = ok.copy(); uniq[ok] &= ~np.isin(idx[ok], list(shared)); uniq[uniq] &= A[idx[uniq]] > 0
j = idx[uniq]; L = lum[uniq]; G = gr[uniq]
rng = np.random.default_rng(1); sub = rng.choice(len(j), 60000, replace=False); js = j[sub]; Ls = L[sub]; Gs = G[sub]
print('testing', len(js), 'records', flush=True)
def rep(y, x, nm):
    m = np.isfinite(x); a, b = np.polyfit(x[m], y[m], 1); r = np.corrcoef(x[m], y[m])[0, 1]; res = y[m] - (a * x[m] + b)
    print('  byte = %8.4f * %-14s + %8.4f   r=%.4f rms=%.2f  lo=%.3f hi=%.3f' % (a, nm, b, r, res.std(), -b / a, (255 - b) / a), flush=True)
rep(Ls, np.log10(L_r[js]), 'logL_r self')
for k in (4, 8, 16, 32, 64, 128):
    dd, ii = tree.query(P[js], k=k)
    S = L_r[ii].sum(1); Sg = L_g[ii].sum(1)
    rep(Ls, np.log10(S), 'logL_r knn%d' % k)
    rep(Gs, -2.5 * np.log10(Sg / S), 'g-r knn%d' % k)
    rep(Gs, (PH[ii, 4] - PH[ii, 5]).mean(1), 'g-r mean knn%d' % k)
for rad in (0.5, 1.0, 2.0, 4.0):
    S = np.zeros(len(js)); Sg = np.zeros(len(js))
    for q, i in enumerate(js):
        ii = tree.query_ball_point(P[i], rad); S[q] = L_r[ii].sum(); Sg[q] = L_g[ii].sum()
    rep(Ls, np.log10(S), 'logL_r r<%.1f' % rad); rep(Gs, -2.5 * np.log10(Sg / S), 'g-r r<%.1f' % rad)
# density-like: log(L_self * n_neighbours)?  and mass-weighted?
dd, ii = tree.query(P[js], k=32); h32 = dd[:, -1]
rep(Ls, np.log10(L_r[js] / h32 ** 3), 'logL_r/h32^3'); rep(Ls, np.log10(1. / h32 ** 3), 'log 1/h32^3')
rep(Ls, np.log10(M[js]), 'logM self'); rep(Ls, A[js], 'age self'); rep(Gs, A[js], 'age self(gr)')
