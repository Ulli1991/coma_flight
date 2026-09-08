# Builds /ptmp/uli/coma_cubes/label_tracers.npz: the 1000 most bound DM particles of the z = 0 subhaloes of
# the named galaxies (NGC 4889 = subhalo 0, NGC 4874 = subhalo 75, NGC 4839 group = subhalo 1, all in
# group 0), from the group catalogue's SubhaloOffsetType.  extract_cubes.py SNAP lab then locates them in
# every snapshot and pack_cubes.py writes the positions into data/ep{S}_meta.json.
#   srun -p interactive -n1 -c4 --mem=24G -t 15 /u/uli/jupyter_env/bin/python3 tools/label_tracers.py
import numpy as np, h5py, glob
OUTD = '/ptmp/uli/coma_cubes'
base = '/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/'
files = sorted(glob.glob(base + 'snapdir_139/snap_139.*.hdf5'), key=lambda p: int(p.split('.')[-2]))
gfiles = sorted(glob.glob(base + 'groups_139/fof_subhalo_tab_139.*.hdf5'), key=lambda p: int(p.split('.')[-2]))
cen = np.load(OUTD + '/cen139.npy'); BOX = 681000.; HALF = 4500.
def wrap(c): d = c - cen; d -= np.round(d / BOX) * BOX; return d
SP = []; SM = []; SLT = []; SGN = []; SOT = []; SST = []
for gf in gfiles:
    with h5py.File(gf) as g:
        if 'Subhalo' in g and 'SubhaloPos' in g['Subhalo']:
            s = g['Subhalo']; SP.append(s['SubhaloPos'][:]); SM.append(s['SubhaloMass'][:]); SLT.append(s['SubhaloLenType'][:])
            SGN.append(s['SubhaloGroupNr'][:]); SOT.append(s['SubhaloOffsetType'][:]); SST.append(s['SubhaloMassType'][:, 4])
SP = np.concatenate(SP).astype(np.float64); SM = np.concatenate(SM); SLT = np.concatenate(SLT); SGN = np.concatenate(SGN); SOT = np.concatenate(SOT); SST = np.concatenate(SST)
npf = []
for fn in files:
    with h5py.File(fn) as f: npf.append(int(f['Header'].attrs['NumPart_ThisFile'][1]))
fstart = np.concatenate([[0], np.cumsum(npf)])
def read_slice(off, cnt, ds):
    out = []
    for i, fn in enumerate(files):
        a, b = max(off, fstart[i]), min(off + cnt, fstart[i + 1])
        if a < b:
            with h5py.File(fn) as f: out.append(f['PartType1/' + ds][a - fstart[i]:b - fstart[i]])
    return np.concatenate(out)
LAB = {'NGC 4889 · BCG': [0, 0, 0], 'NGC 4874': [0.029, -0.071, -0.045], 'NGC 4839 GROUP': [-0.128, 0.517, -0.050]}
dS = wrap(SP); tr = {}
for nm, p in LAB.items():
    target = np.array(p) * HALF; dd = np.linalg.norm(dS - target, axis=1); cand = np.where(dd < 150.)[0]
    order = cand[np.argsort(SM[cand])[::-1]][:4]
    print('label %-16s target %s: candidates within 150 kpc/h (mass, stellar mass, dist):' % (nm, np.round(target)), ', '.join('%d: %.1e %.1e %.0f' % (k, SM[k] * 1e10, SST[k] * 1e10, dd[k]) for k in order))
    k = order[0]; n = int(min(1000, SLT[k, 1]))
    ids = read_slice(SOT[k, 1], n, 'ParticleIDs'); pos = wrap(read_slice(SOT[k, 1], n, 'Coordinates').astype(np.float64))
    med = np.median(pos, axis=0)
    print('   -> subhalo %d (group %d) mass %.2e: %d most bound DM, median at %s, %.1f kpc/h from SubhaloPos, 90%% within %.0f kpc/h' % (
        k, SGN[k], SM[k] * 1e10, n, np.round(med), np.linalg.norm(med - dS[k]), np.percentile(np.linalg.norm(pos - med, axis=1), 90)))
    tr[nm] = ids
np.savez(OUTD + '/label_tracers.npz', names=np.array(list(tr)), **{'ids%d' % i: v for i, v in enumerate(tr.values())})
print('wrote', OUTD + '/label_tracers.npz')
