#!/usr/bin/env python
# Fine time sampling for the Coma Virtual Observatory's time slider: every snapshot 27 .. 139 (z = 1.95 .. 0,
# ~90 Myr apart) as a 128^3 cube of gas mass, mass-weighted log T and r-band star light around the main
# progenitor, plus the tracked positions of the named galaxies.  Run as a job array (run_movie.sbatch):
#   python extract_movie.py SNAP   -> /ptmp/uli/coma_cubes/movie/mv_SNAP.h5
# Centre as in extract_cubes.py: median of the z=0 BCG's most bound DM particles, snapped to the most massive
# FOF group within 500 kpc/h (139: cen139.npy).  pack_cubes.py movie turns these into data/mv/*.u8.gz.
import sys, os, glob, time, numpy as np, h5py
from multiprocessing import Pool
_argv = list(sys.argv); sys.argv = [sys.argv[0], sys.argv[1]]   # extract_cubes reads SNAP from argv at import
import extract_cubes as X

SNAP = X.SNAP; NM = 128; NV = 48; HALF = X.HALF        # NV: the coarse velocity grid for the advected time interpolation
MORE = len(_argv) > 2 and _argv[2] == 'more'   # `extract_movie.py SNAP more`: add gas_mv* (48^3) and the galaxy table to an existing mv file
GALONLY = len(_argv) > 2 and _argv[2] == 'gal'  # `extract_movie.py SNAP gal`: only the galaxy table (catalogue read, seconds)
OUT = X.OUT + '/movie'; RAWF = OUT + '/mv_%03d.h5' % SNAP

def do_movie(args):
    fn, cen = args
    out = {}
    with h5py.File(fn, 'r') as f:
        g = f['PartType0']
        d = X.wrap(g['Coordinates'][:].astype(np.float64), cen)
        m = np.max(np.abs(d), axis=1) < HALF
        if m.any():
            d = d[m]; mass = g['Masses'][:][m].astype(np.float64); T, sfr = X.gas_T(g, m)
            lT = np.log10(np.clip(T, 1e2, None))
            out.update(X.grids(d, {'gas_m': mass, 'gas_mlt': mass * lT}, HALF, NM))
        if 'PartType4' in f and 'Coordinates' in f['PartType4']:
            s = f['PartType4']
            d = X.wrap(s['Coordinates'][:].astype(np.float64), cen); age = s['GFM_StellarFormationTime'][:]
            m = (np.max(np.abs(d), axis=1) < HALF) & (age > 0)
            if m.any():
                L = 10 ** (-0.4 * s['GFM_StellarPhotometrics'][:, 5][m].astype(np.float64))
                out.update(X.grids(d[m], {'st_L': L}, HALF, NM))
    return out

def do_vel(args):
    """mass-weighted gas velocity [km/s physical] on the coarse NV^3 grid (all gas), for the advected interpolation"""
    fn, cen = args
    with h5py.File(fn, 'r') as f:
        g = f['PartType0']
        d = X.wrap(g['Coordinates'][:].astype(np.float64), cen)
        m = np.max(np.abs(d), axis=1) < HALF
        if not m.any(): return {}
        d = d[m]; mass = g['Masses'][:][m].astype(np.float64); v = g['Velocities'][:][m].astype(np.float64) * np.sqrt(1. / (1. + X.Z))
    return X.grids(d, {'gas_mV': mass, 'gas_mvx': mass * v[:, 0], 'gas_mvy': mass * v[:, 1], 'gas_mvz': mass * v[:, 2]}, HALF, NV)

def galaxies(cen):
    """the subhaloes inside the cube with M* >= 1e9 Msun: most-bound particle id (to match neighbouring snapshots),
    position [cube units], log10 M* [Msun], g-r, SFR flag, log10 M_BH [Msun] (h = 0.681), velocity [km/s peculiar, 3]"""
    H = 0.681; P = []; MT = []; SF = []; PH = []; BH = []; ID = []; V = []
    for gf in sorted(glob.glob(X.BASE + '/groups_%03d/fof_subhalo_tab_%03d.*.hdf5' % (SNAP, SNAP))):
        with h5py.File(gf, 'r') as g:
            if 'Subhalo' not in g or 'SubhaloPos' not in g['Subhalo']: continue
            S = g['Subhalo']; P.append(S['SubhaloPos'][:]); MT.append(S['SubhaloMassType'][:]); SF.append(S['SubhaloSFR'][:])
            PH.append(S['SubhaloStellarPhotometrics'][:]); BH.append(S['SubhaloBHMass'][:]); ID.append(S['SubhaloIDMostbound'][:]); V.append(S['SubhaloVel'][:])
    if not P: return np.zeros((0, 11))
    P = np.concatenate(P).astype(np.float64); MT = np.concatenate(MT); SF = np.concatenate(SF); PH = np.concatenate(PH); BH = np.concatenate(BH); ID = np.concatenate(ID); V = np.concatenate(V).astype(np.float64)
    d = X.wrap(P, cen); ms = MT[:, 4] * 1e10 / H
    keep = (np.max(np.abs(d), axis=1) < HALF) & (ms >= 1e9)
    out = np.column_stack([ID[keep].astype(np.float64), d[keep] / HALF, np.log10(ms[keep]), PH[keep, 4] - PH[keep, 5], (SF[keep] > 0).astype(float),
                           np.log10(np.maximum(BH[keep] * 1e10 / H, 1.)), V[keep]])
    return out

def centre(cen0):
    """as extract_cubes.py: tracer median, tightened, snapped to the most massive group within 500 kpc/h"""
    info = {}
    if SNAP == 139: return cen0.copy(), info
    ids = np.load(X.OUT + '/tracers139.npy')
    with Pool(X.NPROC) as p: tr = np.concatenate(p.map(X.find_tracers, [(fn, ids) for fn in X.files]))
    dtr = X.wrap(tr, cen0); med = np.median(dtr, axis=0)
    near = np.linalg.norm(dtr - med, axis=1) < 300.
    if near.sum() > 50: med = np.median(dtr[near], axis=0)
    cen_tr = cen0 + med; info['n_tracers'] = int(len(tr)); info['n_tracers_near'] = int(near.sum())
    GP = []; GM = []; GR = []
    for gf in sorted(glob.glob(X.BASE + '/groups_%03d/fof_subhalo_tab_%03d.*.hdf5' % (SNAP, SNAP))):
        with h5py.File(gf, 'r') as g:
            if 'Group' in g and 'GroupPos' in g['Group']:
                GP.append(g['Group/GroupPos'][:]); GM.append(g['Group/Group_M_Crit200'][:]); GR.append(g['Group/Group_R_Crit200'][:])
    GP = np.concatenate(GP).astype(np.float64); GM = np.concatenate(GM); GR = np.concatenate(GR)
    dg = np.linalg.norm(X.wrap(GP, cen_tr), axis=1); cand = np.where(dg < 500.)[0]
    if len(cand):
        j = cand[np.argmax(GM[cand])]; cen = GP[j]
        info['group_M200'] = float(GM[j]); info['group_R200'] = float(GR[j]); info['group_dist'] = float(dg[j])
    else:
        cen = cen_tr; info['group_M200'] = 0.0; info['group_R200'] = 0.0; info['group_dist'] = -1.0
    return cen, info

if __name__ == '__main__':
    t0 = time.time(); os.makedirs(OUT, exist_ok=True)
    if GALONLY:
        with h5py.File(RAWF, 'r') as f: cen = f.attrs['cen'][:]
        gal = galaxies(cen)
        with h5py.File(RAWF, 'a') as o:
            if 'gal' in o: del o['gal']
            o.create_dataset('gal', data=gal)
        print('galaxy table of', RAWF, len(gal), '%.0f s' % (time.time() - t0), flush=True); sys.exit(0)
    if MORE:   # add the coarse velocity grid and the galaxy table to an existing file (its centre is reused)
        with h5py.File(RAWF, 'r') as f: cen = f.attrs['cen'][:]
        acc = {}
        with Pool(X.NPROC) as p:
            for part in p.imap_unordered(do_vel, [(fn, cen) for fn in X.files]): X.add(acc, part)
        gal = galaxies(cen)
        with h5py.File(RAWF, 'a') as o:
            for k, v in acc.items():
                if k in o: del o[k]
                o.create_dataset(k, data=v, compression='gzip', compression_opts=1)
            if 'gal' in o: del o['gal']
            o.create_dataset('gal', data=gal)
        print('added to', RAWF, sorted(acc), 'galaxies', len(gal), '%.0f s' % (time.time() - t0), flush=True); sys.exit(0)
    cen0 = np.load(X.OUT + '/cen139.npy'); cen, info = centre(cen0)
    print('snap', SNAP, 'z=%.3f' % X.Z, 'centre', cen, 'offset', X.wrap(cen[None], cen0)[0], info, flush=True)
    acc = {}
    with Pool(X.NPROC) as p:
        for part in p.imap_unordered(do_movie, [(fn, cen) for fn in X.files]): X.add(acc, part)
    lab, names = X.track_labels(cen) if os.path.exists(X.OUT + '/label_tracers.npz') else ({}, None)
    with h5py.File(RAWF, 'w') as o:
        o.attrs['snap'] = SNAP; o.attrs['z'] = X.Z; o.attrs['cen'] = cen; o.attrs['cen139'] = cen0; o.attrs['half'] = HALF; o.attrs['n'] = NM
        for k, v in info.items(): o.attrs[k] = v
        for k, v in acc.items(): o.create_dataset(k, data=v, compression='gzip', compression_opts=1)
        for k, v in lab.items(): o.create_dataset(k, data=v)
        if names: o.attrs['lab_names'] = np.array([s.encode('utf-8') for s in names])
    print('wrote', RAWF, sorted(acc), '%.0f s' % (time.time() - t0), flush=True)
