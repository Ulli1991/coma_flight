# Where the Arepo shock-finder cells are at z = 0, by Mach number, temperature and radius.  Result (Sept
# 2026): no shocks in star-forming cells; 42k cells with Mach >= 10 lie beyond a third of the cube at
# 1e5-1e7 K (the accretion shock) and 22 inside; Mach 2-5 cells inside are the ICM merger shocks at
# 1e7-6e7 K.  Justifies dropping the shader's temperature mask.
#   srun -p interactive -n1 -c4 --mem=40G -t 15 /u/uli/jupyter_env/bin/python3 tools/shock_cells.py
import numpy as np, h5py, glob
base = '/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/snapdir_139/'
files = sorted(glob.glob(base + 'snap_139.*.hdf5'), key=lambda p: int(p.split('.')[-2]))
cen = np.load('/ptmp/uli/coma_cubes/cen139.npy'); BOX = 681000.; HALF = 4500.
KB = 1.3807e-16; MP = 1.6726e-24; GAMMA = 5. / 3.
M = []; T = []; R = []; S = []; E = []; RHO = []
for fn in files:
    with h5py.File(fn) as f:
        g = f['PartType0']; c = g['Coordinates'][:].astype(np.float64); d = c - cen; d -= np.round(d / BOX) * BOX
        m = np.max(np.abs(d), axis=1) < HALF
        mach = g['Machnumber'][:][m]; k = mach > 1.
        if not k.any(): continue
        u = g['InternalEnergy'][:][m][k]; xe = g['ElectronAbundance'][:][m][k]; XH = g['GFM_Metals'][:, 0][m][k]
        mu = 4. / (1. + 3. * XH + 4. * XH * xe); T.append((GAMMA - 1.) * u * 1e10 * mu * MP / KB)
        M.append(mach[k]); R.append(np.linalg.norm(d[m][k], axis=1) / HALF); S.append(g['StarFormationRate'][:][m][k])
        E.append(g['EnergyDissipation'][:][m][k]); RHO.append(g['Density'][:][m][k])
M = np.concatenate(M); T = np.concatenate(T); R = np.concatenate(R); S = np.concatenate(S); E = np.concatenate(E); RHO = np.concatenate(RHO)
print('cells with Mach>1 in the cube:', len(M), ' of which SF:', (S > 0).sum(), ' Mach>=100 SF:', ((M >= 100) & (S > 0)).sum())
print('Mach   n cells  | T<1e5   1e5-1e6  1e6-1e7  >1e7  | r<0.33  0.33-0.6  >0.6 (r in half=4500) | median rho, ed')
for lo, hi in ((1, 1.5), (1.5, 2), (2, 3), (3, 5), (5, 10), (10, 20), (20, 100), (100, 1e9)):
    k = (M >= lo) & (M < hi); n = k.sum()
    tb = [((T[k] < 1e5)).sum(), ((T[k] >= 1e5) & (T[k] < 1e6)).sum(), ((T[k] >= 1e6) & (T[k] < 1e7)).sum(), (T[k] >= 1e7).sum()]
    rb = [(R[k] < .33).sum(), ((R[k] >= .33) & (R[k] < .6)).sum(), (R[k] >= .6).sum()]
    print('%4g-%-4g %8d | %7d %7d %7d %7d | %7d %7d %7d | %.2e %.2e' % (lo, hi, n, *tb, *rb, np.median(RHO[k]) if n else 0, np.median(E[k]) if n else 0))
# the accretion shock candidates: Mach>=10 outside 0.33 -- their temperature
k = (M >= 10) & (R >= .33); print('Mach>=10 outside r=0.33: %d cells, T pct 10/50/90: %s' % (k.sum(), np.percentile(T[k], [10, 50, 90]).round(0)))
k = (M >= 10) & (R < .33); print('Mach>=10 inside  r=0.33: %d cells, T pct 10/50/90: %s, rho pct 50/90: %s' % (k.sum(), np.percentile(T[k], [10, 50, 90]).round(0), np.percentile(RHO[k], [50, 90])))
k = (M >= 2) & (M < 5) & (R < .33); print('Mach 2-5 inside r=0.33: %d cells, T pct 10/50/90: %s' % (k.sum(), np.percentile(T[k], [10, 50, 90]).round(0)))
