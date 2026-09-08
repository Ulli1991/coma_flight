#!/usr/bin/env python
# One-pass cube extraction for the Coma Virtual Observatory (run on raven via run_cubes.sbatch).
#   python extract_cubes.py SNAP            -> /ptmp/uli/coma_cubes/raw_SNAP.h5
#   python extract_cubes.py SNAP sp+hires   -> add only the named parts to an existing raw_SNAP.h5
#                                              (reuses its centre; parts: gas dm stars sp hires)
# Raw float32 cubes (192^3 over +-4500 ckpc/h around the main progenitor of the z=0 BCG;
# snapshot 139 additionally 384^3 + 192^3 inner +-1200 ckpc/h; snapshots 109, 121, 139 the
# 384^3 gas mass and star light, "hires") for gas mass, temperature, X-ray proxies, Arepo
# shock-finder Mach number (max in voxel), shock energy dissipation, dark matter mass, r-band
# star light and g-r colour, plus the star particles brighter than SP_MINLOGL ("sp": position,
# log L_r, g-r, formation time) for the galaxy sprites.  pack_cubes.py turns them into the
# u8.gz / bin.gz files the page loads.  Centre for SNAP != 139: median position of the 3000 most
# bound DM particles of the z=0 BCG (tracers139.npy), snapped to the nearest massive group.
import sys, os, glob, time, numpy as np, h5py
from multiprocessing import Pool

SNAP = int(sys.argv[1])
PARTS = sys.argv[2].replace(',', '+').split('+') if len(sys.argv) > 2 else ['gas', 'dm', 'stars', 'sp']
BASE = '/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output'
OUT = '/ptmp/uli/coma_cubes'
HALF = 4500.0; N = 192; N2 = 384; IH = 1200.0; NI = 192
FULL = SNAP == 139                                   # inner cubes + 384^3 shocks / DM
HIRES = SNAP in (109, 121, 139) or 'hires' in PARTS  # 384^3 gas mass + star light (the late epochs)
SP_MINLOGL = 3.0                                     # star particles kept for the sprites (log10 L_r, mag zero point)
KB = 1.3807e-16; MP = 1.6726e-24; GAMMA = 5. / 3.
NPROC = int(os.environ.get('SLURM_CPUS_PER_TASK', '16'))
files = sorted(glob.glob(BASE + '/snapdir_%03d/snap_%03d.*.hdf5' % (SNAP, SNAP)), key=lambda p: int(p.split('.')[-2]))
if os.environ.get('CUBES_FILES'): files = [files[int(i)] for i in os.environ['CUBES_FILES'].split(',')]   # quick test subset
with h5py.File(files[0], 'r') as f:
    BOX = float(f['Header'].attrs['BoxSize']); Z = float(f['Header'].attrs['Redshift']); MT = f['Header'].attrs['MassTable'][:]

def wrap(c, cen):
    d = c - cen; d -= np.round(d / BOX) * BOX; return d

def flat(d, half, n):
    i = np.clip(((d + half) / (2 * half) * n).astype(np.int64), 0, n - 1)
    return (i[:, 0] * n + i[:, 1]) * n + i[:, 2]

def hist(idx, w, n):
    return np.bincount(idx, w, minlength=n ** 3).astype(np.float32).reshape(n, n, n)

def hmax(idx, v, n):
    g = np.zeros(n ** 3, np.float32); np.maximum.at(g, idx, v.astype(np.float32)); return g.reshape(n, n, n)

def grids(d, w, half, n, extra=None):
    """mass-like histograms for a dict of weights, plus max-in-cell fields in `extra`"""
    idx = flat(d, half, n); out = {}
    for k, v in w.items(): out[k] = hist(idx, v, n)
    for k, v in (extra or {}).items(): out[k] = hmax(idx, v, n)
    return out

def add(acc, part):
    for k, v in part.items():
        if k not in acc: acc[k] = v
        elif k.startswith('max_'): np.maximum(acc[k], v, out=acc[k])
        elif k.startswith('sp_'): acc[k] = np.concatenate([acc[k], v])
        else: acc[k] += v

# ---------------------------------------------------------------- tracers -> centre
def find_tracers(args):
    fn, ids = args
    with h5py.File(fn, 'r') as f:
        pid = f['PartType1/ParticleIDs'][:]
        m = np.isin(pid, ids)
        return f['PartType1/Coordinates'][:][m] if m.any() else np.zeros((0, 3))

# ---------------------------------------------------------------- per-file workers
def do_gas(args):
    fn, cen = args
    with h5py.File(fn, 'r') as f:
        g = f['PartType0']
        d = wrap(g['Coordinates'][:].astype(np.float64), cen)
        m = np.max(np.abs(d), axis=1) < HALF
        if not m.any(): return {}
        d = d[m]
        mass = g['Masses'][:][m].astype(np.float64); rho = g['Density'][:][m].astype(np.float64)
        u = g['InternalEnergy'][:][m]; xe = g['ElectronAbundance'][:][m]; XH = g['GFM_Metals'][:, 0][m]
        sfr = g['StarFormationRate'][:][m]; mach = g['Machnumber'][:][m]; ed = g['EnergyDissipation'][:][m].astype(np.float64)
    mu = 4. / (1. + 3. * XH + 4. * XH * xe)
    T = (GAMMA - 1.) * u * 1e10 * mu * MP / KB
    T[sfr > 0] = 1e4                                        # effective-EOS cells: treat as cold ISM
    lT = np.log10(np.clip(T, 1e2, None))
    hot = T > 1e6
    ism = sfr > 0                                           # keep the effective-EOS ISM out of the shock fields (the
    mach = np.where(ism, 0., mach); ed = ed * ~ism           # shock finder already skips it: a no-op guard at z=0)
    w = {'gas_m': mass, 'gas_mlt': mass * lT,
         'gas_xr': mass * rho * np.sqrt(T) * hot,           # bremsstrahlung: rho^2 V sqrt(T)
         'gas_em': mass * rho * hot,                        # emission measure, flat cooling function
         'gas_ed': ed}
    out = grids(d, w, HALF, N, {'max_mach': mach})
    if HIRES: out.update(grids(d, {'gas_m384': mass}, HALF, N2))
    if FULL:
        o2 = grids(d, {'gas_ed384': ed}, HALF, N2, {'max_mach384': mach}); out.update(o2)
        mi = np.max(np.abs(d), axis=1) < IH
        o3 = grids(d[mi], {'gas_mI': mass[mi], 'gas_edI': ed[mi]}, IH, NI, {'max_machI': mach[mi]}); out.update(o3)
    return out

def do_dm(args):
    fn, cen = args
    out = {}
    with h5py.File(fn, 'r') as f:
        for pt in (1, 2):
            k = 'PartType%d' % pt
            if k not in f or 'Coordinates' not in f[k]: continue
            d = wrap(f[k + '/Coordinates'][:].astype(np.float64), cen)
            m = np.max(np.abs(d), axis=1) < HALF
            if not m.any(): continue
            d = d[m]
            wgt = f[k + '/Masses'][:][m].astype(np.float64) if 'Masses' in f[k] else np.full(len(d), float(MT[pt]))
            add(out, grids(d, {'dm_m': wgt}, HALF, N))
            if FULL:
                add(out, grids(d, {'dm_m384': wgt}, HALF, N2))
                mi = np.max(np.abs(d), axis=1) < IH
                add(out, grids(d[mi], {'dm_mI': wgt[mi]}, IH, NI))
    return out

def do_stars(args):
    fn, cen = args
    with h5py.File(fn, 'r') as f:
        if 'PartType4' not in f or 'Coordinates' not in f['PartType4']: return {}
        s = f['PartType4']
        d = wrap(s['Coordinates'][:].astype(np.float64), cen)
        age = s['GFM_StellarFormationTime'][:]
        m = (np.max(np.abs(d), axis=1) < HALF) & (age > 0)          # drop wind-phase particles
        if not m.any(): return {}
        d = d[m]; ph = s['GFM_StellarPhotometrics'][:][m]            # U B V K g r i z absolute mags
    L = 10 ** (-0.4 * ph[:, 5].astype(np.float64))                   # r-band luminosity (arbitrary zero point)
    gr = (ph[:, 4] - ph[:, 5]).astype(np.float64)
    out = grids(d, {'st_L': L, 'st_Lc': L * gr}, HALF, N)
    if HIRES: out.update(grids(d, {'st_L384': L, 'st_Lc384': L * gr}, HALF, N2))
    if FULL:
        mi = np.max(np.abs(d), axis=1) < IH
        out.update(grids(d[mi], {'st_LI': L[mi], 'st_LcI': L[mi] * gr[mi]}, IH, NI))
    return out

def do_sp(args):
    """star particles for the galaxy sprites: position (ckpc/h from the centre), log10 L_r, g-r, formation a"""
    fn, cen = args
    with h5py.File(fn, 'r') as f:
        if 'PartType4' not in f or 'Coordinates' not in f['PartType4']: return {}
        s = f['PartType4']
        d = wrap(s['Coordinates'][:].astype(np.float64), cen)
        age = s['GFM_StellarFormationTime'][:]
        m = (np.max(np.abs(d), axis=1) < HALF) & (age > 0)
        if not m.any(): return {}
        ph = s['GFM_StellarPhotometrics'][:][m]
    lr = -0.4 * ph[:, 5]; keep = lr >= SP_MINLOGL
    return {'sp_pos': d[m][keep].astype(np.float32), 'sp_lr': lr[keep].astype(np.float32),
            'sp_gr': (ph[:, 4] - ph[:, 5])[keep].astype(np.float32), 'sp_age': age[m][keep].astype(np.float32)}

if __name__ == '__main__':
    t0 = time.time(); os.makedirs(OUT, exist_ok=True)
    cen0 = np.load(OUT + '/cen139.npy')
    info = {}; RAWF = OUT + '/raw_%03d.h5' % SNAP
    APPEND = len(sys.argv) > 2 and os.path.exists(RAWF)
    if APPEND:
        with h5py.File(RAWF, 'r') as f: cen = f.attrs['cen'][:]; Z = float(f.attrs['z'])
        print('snap', SNAP, 'adding', PARTS, 'to', RAWF, 'centre', cen, flush=True)
    elif FULL:
        cen = cen0.copy()
    else:
        ids = np.load(OUT + '/tracers139.npy')
        with Pool(NPROC) as p: tr = np.concatenate(p.map(find_tracers, [(fn, ids) for fn in files]))
        dtr = wrap(tr, cen0)
        med = np.median(dtr, axis=0); cen_tr = cen0 + med
        # tighten: median of the tracers within 300 kpc/h of the first estimate
        near = np.linalg.norm(dtr - med, axis=1) < 300.
        if near.sum() > 50: med = np.median(dtr[near], axis=0); cen_tr = cen0 + med
        info['n_tracers'] = int(len(tr)); info['n_tracers_near'] = int(near.sum()); info['cen_tracer'] = cen_tr
        # snap to the most massive FOF group within 500 kpc/h of the tracer centre
        GP = []; GM = []; GR = []
        for gf in sorted(glob.glob(BASE + '/groups_%03d/fof_subhalo_tab_%03d.*.hdf5' % (SNAP, SNAP))):
            with h5py.File(gf, 'r') as g:
                if 'Group' in g and 'GroupPos' in g['Group']:
                    GP.append(g['Group/GroupPos'][:]); GM.append(g['Group/Group_M_Crit200'][:]); GR.append(g['Group/Group_R_Crit200'][:])
        GP = np.concatenate(GP).astype(np.float64); GM = np.concatenate(GM); GR = np.concatenate(GR)
        dg = np.linalg.norm(wrap(GP, cen_tr), axis=1); cand = np.where(dg < 500.)[0]
        if len(cand):
            j = cand[np.argmax(GM[cand])]; cen = GP[j]
            info['group_M200'] = float(GM[j]); info['group_R200'] = float(GR[j]); info['group_dist'] = float(dg[j])
        else:
            cen = cen_tr; info['group_M200'] = 0.0; info['group_R200'] = 0.0; info['group_dist'] = -1.0
        print('snap', SNAP, 'z=%.3f' % Z, 'tracers', len(tr), 'centre', cen, 'offset from z=0 centre', wrap(cen[None], cen0)[0], info, flush=True)
    acc = {}
    todo = [(n, f) for n, f in (('gas', do_gas), ('dm', do_dm), ('stars', do_stars), ('sp', do_sp))
            if n in PARTS or ('hires' in PARTS and n in ('gas', 'stars'))]
    for name, fun in todo:
        with Pool(NPROC) as p:
            for part in p.imap_unordered(fun, [(fn, cen) for fn in files]): add(acc, part)
        print(name, 'done', '%.0f s' % (time.time() - t0), flush=True)
    if 'sp_lr' in acc: print('sp: %d star particles with log L_r >= %.1f' % (len(acc['sp_lr']), SP_MINLOGL), flush=True)
    with h5py.File(RAWF, 'a' if APPEND else 'w') as o:
        if not APPEND:
            o.attrs['snap'] = SNAP; o.attrs['z'] = Z; o.attrs['cen'] = cen; o.attrs['cen139'] = cen0; o.attrs['half'] = HALF
            o.attrs['inner_half'] = IH; o.attrs['box'] = BOX
            for k, v in info.items(): o.attrs[k] = v
        for k, v in acc.items():
            if k in o: del o[k]
            o.create_dataset(k, data=v, compression='gzip', compression_opts=1)
    print('wrote', RAWF, 'keys', sorted(acc), '%.0f s' % (time.time() - t0), flush=True)
