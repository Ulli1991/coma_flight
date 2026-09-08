# How far the shipped galaxy sprites (data/stars.bin.gz) can be reproduced from snapshot 139: measures the
# centre offset of the records against the star particles (0.05 kpc/h, i.e. none), matches records to
# particles within 0.3 kpc/h and fits the lum / g-r bytes against the particle's own photometry.
# Result (Sept 2026): r = 0.70 (lum vs log L_r) and 0.79 (g-r), strongly non-linear -> pack_cubes.py
# calib_sp matches distributions instead (see fit_star_sprites_neighbours.py for the neighbourhood test).
#   srun -p interactive -n1 -c8 --mem=60G -t 20 /u/uli/jupyter_env/bin/python3 tools/fit_star_sprites.py
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
t0 = time.time()
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
print('star particles', len(P), '%.0f s' % (time.time() - t0), flush=True)
tree = cKDTree(P)
# offset: displacement histogram between bright shipped records and all particles within 6 kpc/h
b = np.argsort(lum)[::-1][:3000]
disp = []
for i in b:
    ii = tree.query_ball_point(xyz[i], 6.0)
    if ii: disp.append(P[ii] - xyz[i])
disp = np.concatenate(disp); print('displacement samples', len(disp))
Hh, edges = np.histogramdd(disp, bins=120, range=[[-6, 6]] * 3)
k = np.unravel_index(np.argmax(Hh), Hh.shape); pk = [(edges[a][k[a]] + edges[a][k[a] + 1]) / 2 for a in range(3)]
print('peak displacement (particle - record)', np.round(pk, 2), 'count', Hh.max(), 'median count', np.median(Hh[Hh > 0]))
# refine: mean of displacements within 0.3 of the peak
near = np.linalg.norm(disp - pk, axis=1) < 0.3; off = disp[near].mean(0)
print('refined offset', np.round(off, 3), 'n', near.sum())
np.save('/ptmp/uli/coma_cubes/star_offset.npy', off)
# rematch with the offset
dist, idx = tree.query(xyz + off, distance_upper_bound=0.3)
ok = np.isfinite(dist)
print('records matched within 0.3:', ok.sum(), 'of', n, 'dist pct', np.percentile(dist[ok], [50, 90, 99]))
# unique (nearest particle not shared)
u, cnt = np.unique(idx[ok], return_counts=True); shared = set(u[cnt > 1])
uniq = ok.copy(); uniq[ok] &= ~np.isin(idx[ok], list(shared))
j = idx[uniq]; ph = PH[j]; age = A[j]; mass = M[j]; L = lum[uniq]; G = gr[uniq]
print('unique', uniq.sum(), 'wind', (age <= 0).sum())
good = age > 0
ph = ph[good]; L = L[good]; G = G[good]; age = age[good]; mass = mass[good]
Lr = -0.4 * ph[:, 5]; g_r = ph[:, 4] - ph[:, 5]
def rep(y, x, nm):
    a, b = np.polyfit(x, y, 1); r = np.corrcoef(x, y)[0, 1]; res = y - (a * x + b)
    print('  byte = %.4f * %s + %.4f   r=%.5f  rms=%.3f  -> lo=%.4f hi=%.4f' % (a, nm, b, r, res.std(), -b / a, (255 - b) / a))
    return a, b
for nm, x in (('logL_r', Lr), ('logL_V', -0.4 * ph[:, 2]), ('logL_g', -0.4 * ph[:, 4]), ('logL_K', -0.4 * ph[:, 7]), ('logM', np.log10(mass))): rep(L, x, nm)
for nm, x in (('g-r', g_r), ('g-i', ph[:, 4] - ph[:, 6]), ('B-V', ph[:, 1] - ph[:, 2]), ('age', age)): rep(G, x, nm)
# quantiles of the unmatched-but-bright question: what is the selection?
allLr = -0.4 * PH[:, 5]; sel = np.zeros(len(P), bool); sel[j[good]] = True
thr = Lr.min(); print('min logL_r among matched %.3f, 1st pct %.3f; particles above min: %d, above 1st pct: %d' % (thr, np.percentile(Lr, 1), (allLr >= thr).sum(), (allLr >= np.percentile(Lr, 1)).sum()))
print('mass: matched min %.3e  pct1 %.3e ; all particles above pct1: %d' % (mass.min(), np.percentile(mass, 1), (M >= np.percentile(mass, 1)).sum()))
# is the selection random? fraction of particles selected vs logL bin
for lo_, hi_ in ((3, 3.5), (3.5, 4), (4, 4.5), (4.5, 5), (5, 5.5), (5.5, 6), (6, 7)):
    mm = (allLr >= lo_) & (allLr < hi_); print('  logL_r %.1f-%.1f: %d particles, %d matched (%.2f%%)' % (lo_, hi_, mm.sum(), (mm & sel).sum(), 100 * (mm & sel).sum() / max(mm.sum(), 1)))
np.savez('/ptmp/uli/coma_cubes/match139.npz', Lr=Lr, g_r=g_r, L=L, G=G, age=age, mass=mass, off=off)
