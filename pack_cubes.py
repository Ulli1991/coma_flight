#!/usr/bin/env python
# Turn the raw cubes of extract_cubes.py into the u8.gz files the page loads.
#   python pack_cubes.py check   -> centre / axis-order check of raw_139.h5 against the shipped rho384
#   python pack_cubes.py calib   -> the check, then fits the u8 scalings of the existing z=0 cubes
#                                   (data/*.u8.gz) against raw_139.h5 and writes cube_scales.json
#   python pack_cubes.py calib_sp -> fits the galaxy-sprite encoding (sp_lum, sp_col, sp_min) of the shipped
#                                   data/stars.bin.gz against the star particles of raw_139.h5 (sp_* datasets)
#   python pack_cubes.py calib_cg -> the cold-gas sprite encoding of data/gas.bin.gz against raw_139.h5 (cg_* datasets)
#   python pack_cubes.py 139     -> data/dm384.u8.gz idm192.u8.gz shock384.u8.gz ishock192.u8.gz
#   python pack_cubes.py kin|met|bvec -> vel192 (RGBA: v vector + log sigma), met192, bvec192 (+ ep*_ versions)
#   python pack_cubes.py outer   -> data/outer192.u8.gz (RG: matter + gas over +-30 Mpc/h, z=0)
#   python pack_cubes.py movie   -> data/mv/mv_SNAP.u8.gz (128^3 RGB) + data/movie.json from extract_movie.py
#   python pack_cubes.py galaxies -> data/galaxies.json (z=0 subhaloes for the galaxy inspector)
#   python pack_cubes.py temp    -> temp192.u8.gz + ep*_temp192.u8.gz on the 3e5 .. 1.1e8 K scale (then re-run movie)
#   python pack_cubes.py mag     -> data/mag192.u8.gz + data/ep*_mag192.u8.gz for every raw file with the magnetic
#                                   field (extract_cubes.py SNAP mag); the log scale is fitted once at z=0 (`mag`)
#   python pack_cubes.py 27      -> data/ep027_{pk,xray,temp,dm,shock}192.u8.gz (same scales as z=0)
#                                   + ep027_stars.bin.gz (galaxy sprites), ep027_gas.bin.gz (cold-gas sprites),
#                                   ep027_meta.json (black holes + tracked label positions) + ep027_pk384.u8.gz
#                                   where the raw file has the 384^3 grids (snapshots 109, 121)
# All cubes are C-ordered (x,y,z) uint8; the page samples them as uv.zyx.
import sys, os, json, gzip, glob, numpy as np, h5py
from scipy.ndimage import gaussian_filter, maximum_filter
RAW = '/ptmp/uli/coma_cubes'; DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
SC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cube_scales.json')
LOG8 = np.log10(8.)          # a 192^3 voxel holds 8x the mass of a 384^3 voxel
DMF = 28                     # DM shader: smoothstep(0.12,1.0,dd) -> u8 < 31 is invisible
HALF = 4500.                 # sprite positions: int16 = pos / HALF * 32767
SP_CAP = 1200000             # at most this many sprites per epoch (~7 MB gzipped)

def rd(name, n):
    return np.frombuffer(gzip.decompress(open(os.path.join(DATA, name), 'rb').read()), np.uint8).reshape(n, n, n)

def wr(name, u8):
    path = os.path.join(DATA, name)
    with open(path, 'wb') as f: f.write(gzip.compress(np.ascontiguousarray(u8).tobytes(), 9)); f.flush(); os.fsync(f.fileno())
    print('  wrote', name, '%.1f MB' % (os.path.getsize(path) / 1e6), flush=True)

def u8(q, lo, hi, floor=0, quant=1):
    """floor: values the shader never shows (below its smoothstep threshold) are zeroed, quant: keep 256/quant
    levels -- both only to keep the gzip small (the 384^3 DM cube goes 25 -> 15 MB with 64 levels)"""
    u = (np.clip((q - lo) / (hi - lo), 0, 1) * 255 + 0.5).astype(np.uint8); u[u < floor] = 0
    return (u // quant) * quant if quant > 1 else u

def fit(x, y, name, ok=None):
    """least-squares y = a x + b on the unsaturated range of populated voxels; (lo, hi) with u8 = 255 (x-lo)/(hi-lo)"""
    m = np.isfinite(x) & (y > 4) & (y < 251)
    if ok is not None: m &= ok
    x = x[m].astype(np.float64); y = y[m].astype(np.float64)
    a, b = np.polyfit(x, y, 1); r = np.corrcoef(x, y)[0, 1]
    lo = -b / a; hi = lo + 255. / a
    print('  %-8s lo=%.4f hi=%.4f  r=%.4f  (n=%d)' % (name, lo, hi, r, m.sum()))
    return [float(lo), float(hi), float(r)]

def logq(q):
    return np.log10(np.clip(q, 1e-30, None))

def logT(f, suf=''):
    m = f['gas_m' + suf][:]; return np.where(m > 0, f['gas_mlt' + suf][:] / np.maximum(m, 1e-30), 0.)

def check(f=None):
    """Centre and axis order of raw_139.h5 against the shipped rho384: FFT cross-correlation of the smoothed,
    u8-like stretched log mass histogram for all six axis permutations.  Expect axes (0, 1, 2) at zero shift
    (corr ~0.9; the voxelwise fit in calib() then reaches r = 0.994).  A shift of s voxels means the centre
    is off by s * 23.4 ckpc/h.  (Verification step of the former extract_raven.py, kept here.)"""
    import itertools
    f = f or h5py.File(RAW + '/raw_139.h5', 'r')
    a = rd('rho384.u8.gz', 384).astype(np.float32); a -= a.mean(); n = a.shape[0]
    m = f['gas_m384'][:]; lg = logq(gaussian_filter(m, 1.0)); lo, hi = np.percentile(lg[m > 0], [5, 99.9])
    x = np.clip((lg - lo) / (hi - lo), 0, 1).astype(np.float32); del m, lg
    fa = np.fft.rfftn(a); best = None
    print('centre / axes (shipped rho384 vs smoothed mass histogram, 384^3):')
    for perm in itertools.permutations(range(3)):
        b = np.transpose(x, perm); b = b - b.mean()
        c = np.fft.irfftn(fa * np.conj(np.fft.rfftn(b)), s=a.shape, axes=(0, 1, 2))
        s = np.array(np.unravel_index(np.argmax(c), c.shape)); s[s > n // 2] -= n
        r = float(c.max() / np.sqrt((a * a).sum() * (b * b).sum()))
        print('  axes %s: shift %s voxels, corr %.3f' % (perm, s, r))
        if best is None or r > best[0]: best = (r, perm, s)
    ok = best[1] == (0, 1, 2) and not np.any(best[2])
    print('  -> %s' % ('centre and axis order confirmed (corr %.3f)' % best[0] if ok else
                      'MISMATCH: best axes %s, shift %s voxels = %s ckpc/h; fix cen139.npy / the transpose before packing' % (best[1], best[2], np.round(best[2] * 9000. / n, 1))))
    return ok

def calib():
    f = h5py.File(RAW + '/raw_139.h5', 'r'); S = json.load(open(SC)) if os.path.exists(SC) else {}   # keep sp_* of calib_sp
    if not check(f): sys.exit('centre / axis check failed')
    print('density (384 mass histogram vs rho384; also with 1-voxel smoothing):')
    m384 = f['gas_m384'][:]; ref = rd('rho384.u8.gz', 384)
    fit(logq(m384), ref, 'raw', m384 > 0)                                  # r~0.37: the page cube is smoothed
    S['rho'] = fit(logq(gaussian_filter(m384, 1.0)), ref, 'smooth1', m384 > 0)   # r=0.994
    print('temperature (mass-weighted log T vs temp192):')
    lt = logT(f); S['temp'] = fit(lt, rd('temp192.u8.gz', 192), 'logT', f['gas_m'][:] > 0)   # r=0.93, 1.3e7..1.1e8 K
    print('X-ray (vs xray192):')
    xr = f['gas_xr'][:]; S['xray'] = fit(logq(xr), rd('xray192.u8.gz', 192), 'rho2sqrtT', xr > 0)   # r=0.95
    em = f['gas_em'][:]; fit(logq(em), rd('xray192.u8.gz', 192), 'rho2', em > 0)                     # r=0.93
    print('star light (r-band luminosity per voxel vs slum384, g-r vs scol384):')
    L = f['st_L384'][:]; S['slum'] = fit(logq(L), rd('slum384.u8.gz', 384), 'logL', L > 0)       # r=0.9997, unsmoothed
    gr = np.where(L > 0, f['st_Lc384'][:] / np.maximum(L, 1e-30), np.nan); S['scol'] = fit(gr, rd('scol384.u8.gz', 384), 'g-r')   # r=0.9994
    print('dark matter (384: smooth 1, percentiles 35/99.95 of log, the recipe of the former extract_dm.py):')
    g = gaussian_filter(f['dm_m384'][:], 1.0); lg = logq(g); lo, hi = np.percentile(lg[g > 0], [35, 99.95]); S['dm'] = [float(lo), float(hi), 1.0]
    print('  dm lo=%.4f hi=%.4f' % (lo, hi))
    print('shocks: max Mach per voxel (384) and log dissipation:')
    M = f['max_mach384'][:]; ed = f['gas_ed384'][:]
    print('  voxels with Mach>1: %d (%.2f%%)  Mach pct 50/90/99/99.9: %s' % ((M > 1).sum(), 100 * (M > 1).mean(), np.round(np.percentile(M[M > 1], [50, 90, 99, 99.9]), 2)))
    le = logq(ed[ed > 0]); lo, hi = np.percentile(le, [40, 99.9]); S['ed'] = [float(lo), float(hi), 1.0]
    print('  ed>0 voxels %d, log ed lo=%.3f hi=%.3f' % ((ed > 0).sum(), lo, hi))
    S['mach'] = [1.0, 5.0, 1.0]
    json.dump(S, open(SC, 'w'), indent=1); print('wrote', SC)

def shock_pack(M, ed, S, sh=0.):
    """RG cube: r = Mach (1..5), g = log dissipation.  Shock surfaces are 1 voxel thick and the ray marcher
    steps ~3 voxels, so they are spread with a 1-voxel gaussian and renormalised: a sheet keeps its peak
    (0.4 x 2.5), isolated single voxels (ISM shocks inside galaxies, Mach 100+) drop to 0.16 of theirs."""
    lo, hi = S['mach'][:2]
    m = np.clip(gaussian_filter(np.clip((M - lo) / (hi - lo), 0, 1), 1.0) * 2.5, 0, 1)
    e = np.clip((logq(ed) - (S['ed'][0] + sh)) / (S['ed'][1] - S['ed'][0]), 0, 1) * (ed > 0); e = np.clip(gaussian_filter(e, 1.0) * 2.5, 0, 1)
    r = (m * 255 + .5).astype(np.uint8); r[r < 10] = 0; r = (r // 4) * 4                # below smoothstep(0.04,..): invisible
    g = (e * 255 + .5).astype(np.uint8); g[r == 0] = 0; g = (g // 16) * 16              # brightness modulation only: 16 levels
    return np.stack([r, g], axis=-1)

def rd_stars(name):
    """the page's sprite record: int16 xyz (/32767*HALF ckpc/h), bytes lum, g-r, h, pad -> (n, 10) uint8"""
    raw = gzip.decompress(open(name, 'rb').read()); n = len(raw) // 10
    return np.frombuffer(raw, np.uint8).reshape(n, 10)

def stars_pack(f, S):
    """galaxy sprites of one snapshot in the stars.bin.gz record format, on the calibrated z=0 encoding:
    star particles with log L_r >= sp_min (the z=0 threshold), brightest SP_CAP if there are more"""
    lr = f['sp_lr'][:]; keep = lr >= S['sp_min']
    if keep.sum() > SP_CAP: keep = lr >= np.sort(lr)[-SP_CAP]
    pos = f['sp_pos'][:][keep]; lr = lr[keep]; gr = f['sp_gr'][:][keep]; n = len(lr)
    rec = np.zeros((n, 10), np.uint8)
    rec[:, :6] = np.clip(np.round(pos / HALF * 32767.), -32768, 32767).astype('<i2').view(np.uint8).reshape(n, 6)
    rec[:, 6] = qmap(lr, S['sp_lum']); rec[:, 7] = qmap(gr, S['sp_col']); rec[:, 8] = 90
    print('  sprites: %d star particles (log L_r >= %.3f)' % (n, lr.min()), flush=True)
    return rec

QK = [0, 0.1, 0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98, 99, 99.5, 99.9, 99.99, 100]   # knot quantiles

def qmap(x, table):
    """monotone piecewise-linear map x -> byte through the (x, byte) knots of cube_scales.json"""
    t = np.asarray(table); return (np.interp(x, t[:, 0], t[:, 1]) + 0.5).astype(np.uint8)

def qq_table(x, y, name):
    """quantile-quantile lookup: the byte distribution of the shipped file as a function of the particle
    quantity, matched quantile by quantile (both sorted); knots at QK, ties in x broken by a tiny slope"""
    x = np.sort(x.astype(np.float64)); y = np.sort(y.astype(np.float64))
    kx = np.percentile(x, QK); ky = np.percentile(y, QK)
    kx = np.maximum.accumulate(kx + np.arange(len(kx)) * 1e-6)
    print('  %-4s knots (x -> byte):' % name, ' '.join('%.3g:%.0f' % (a, b) for a, b in zip(kx, ky)))
    return [[float(a), float(b)] for a, b in zip(kx, ky)]

def calib_sp():
    """The shipped stars.bin.gz (built on viper) holds the 939637 brightest star particles; its lum and g-r
    bytes are noisy, non-linear functions of the particle's own r-band luminosity and g-r (r ~ 0.7, 0.8 per
    particle).  Match the *distributions* instead: luminosity threshold = the n-th brightest particle at
    z=0, and a quantile-quantile lookup table byte(log L_r), byte(g-r) over those particles, so an epoch
    file drawn with the same shader looks like the z=0 file at z=0 (same byte histograms)."""
    S = json.load(open(SC)); f = h5py.File(RAW + '/raw_139.h5', 'r')
    ship = rd_stars(os.path.join(DATA, 'stars.bin.gz')); n = len(ship)
    lr = f['sp_lr'][:]; gr = f['sp_gr'][:]
    order = np.argsort(lr)[::-1][:n]; thr = float(lr[order[-1]])
    print('shipped sprites %d, threshold log L_r = %.4f (of %d particles with log L_r >= %.1f)' % (n, thr, len(lr), lr.min()))
    S['sp_lum'] = qq_table(lr[order], ship[:, 6], 'lum'); S['sp_col'] = qq_table(gr[order], ship[:, 7], 'g-r'); S['sp_min'] = thr; S['sp_n'] = n
    json.dump(S, open(SC, 'w'), indent=1); print('wrote', SC)
    # how close is the regenerated z=0 file to the shipped one?  histogram of the bytes and a 96^3 map of the
    # sprite light (0.08 + 2.2 (byte/255)^2.5, as the shader) -- correlation of the two maps
    new = stars_pack(f, S); open(RAW + '/stars139_regen.bin.gz', 'wb').write(gzip.compress(new.tobytes(), 9))
    for k, nm in ((6, 'lum'), (7, 'g-r')):
        print('  %s byte pct 1/10/50/90/99: shipped %s  regenerated %s' % (nm, np.percentile(ship[:, k], [1, 10, 50, 90, 99]), np.percentile(new[:, k], [1, 10, 50, 90, 99])))
    def lmap(rec):
        xyz = rec[:, :6].copy().view('<i2').reshape(len(rec), 3).astype(np.float64) / 32767.
        i = np.clip(((xyz + 1) / 2 * 96).astype(np.int64), 0, 95); idx = (i[:, 0] * 96 + i[:, 1]) * 96 + i[:, 2]
        return np.bincount(idx, 0.08 + 2.2 * (rec[:, 6] / 255.) ** 2.5, minlength=96 ** 3)
    a, b = lmap(ship), lmap(new); m = (a > 0) | (b > 0)
    print('  96^3 sprite-light maps: corr(log) = %.4f over %d cells, total light ratio new/shipped = %.3f' % (np.corrcoef(np.log10(a[m] + .1), np.log10(b[m] + .1))[0, 1], m.sum(), b.sum() / a.sum()))
    print('  wrote', RAW + '/stars139_regen.bin.gz', '(not shipped; copy over data/stars.bin.gz to use the regenerated z=0 sprites)')

BH_LO, BH_HI, BH_MIN = 5.833, 10.252, 0.01583   # D_BH mass byte: m = (log10 M_BH[Msun/h] - lo)/(hi - lo), floor at the seed mass (fit r = 1.000)

def calib_cg():
    """cold-gas sprite encoding of the shipped gas.bin.gz (8-byte records: int16 xyz, byte log rho, byte SF flag)
    against the cold cells of raw_139.h5: the byte is linear in log10 rho (r = 0.98, lo -3.21 hi -0.88 with a floor
    of 9); a quantile table reproduces its histogram exactly.  Shipped set = the 200k densest cold / SF cells."""
    S = json.load(open(SC)); f = h5py.File(RAW + '/raw_139.h5', 'r')
    raw = gzip.decompress(open(os.path.join(DATA, 'gas.bin.gz'), 'rb').read()); n = len(raw) // 8
    ship = np.frombuffer(raw, np.uint8).reshape(n, 8)
    lr = f['cg_lrho'][:]; order = np.argsort(lr)[::-1][:n]
    print('shipped cold-gas sprites %d, raw cold cells %d, threshold log rho = %.3f' % (n, len(lr), lr[order[-1]]))
    S['cg_lrho'] = qq_table(lr[order], ship[:, 6], 'lrho'); S['cg_n'] = n
    json.dump(S, open(SC, 'w'), indent=1); print('wrote', SC)

def gas_pack(f, S):
    lr = f['cg_lrho'][:]; keep = np.argsort(lr)[::-1][:S['cg_n']]
    pos = f['cg_pos'][:][keep]; n = len(keep); rec = np.zeros((n, 8), np.uint8)
    rec[:, :6] = np.clip(np.round(pos / HALF * 32767.), -32768, 32767).astype('<i2').view(np.uint8).reshape(n, 6)
    rec[:, 6] = qmap(lr[keep], S['cg_lrho']); rec[:, 7] = f['cg_sf'][:][keep].astype(np.uint8) * 255
    print('  cold gas: %d cells (log rho >= %.2f, SF %d)' % (n, lr[keep].min(), int(rec[:, 7].sum() // 255)), flush=True)
    return rec

def meta_json(f, name):
    """black holes (x, y, z, m as in the page's D_BH) and the tracked label positions of one snapshot"""
    out = {}
    if 'bh_mass' in f:
        m = np.clip((np.log10(f['bh_mass'][:].astype(np.float64) * 1e10) - BH_LO) / (BH_HI - BH_LO), BH_MIN, 1.)
        p = f['bh_pos'][:] / HALF; out['bh'] = [round(float(v), 5) for row in np.column_stack([p, m]) for v in row]
    if 'lab_pos' in f:
        names = [s.decode() for s in f.attrs['lab_names']]; lp = f['lab_pos'][:] / HALF; ln = f['lab_n'][:]
        out['labels'] = [{'n': nm, 'p': [round(float(v), 4) for v in lp[i]], 'n_tracers': int(ln[i])} for i, nm in enumerate(names) if np.isfinite(lp[i]).all()]
    json.dump(out, open(os.path.join(DATA, name), 'w'), separators=(',', ':'))
    print('  wrote %s: %d black holes, labels %s' % (name, len(out.get('bh', [])) // 4, [(l['n'], l['p'], l['n_tracers']) for l in out.get('labels', [])]), flush=True)

def pack139():
    S = json.load(open(SC)); f = h5py.File(RAW + '/raw_139.h5', 'r')
    lo, hi = S['dm'][:2]
    wr('dm384.u8.gz', u8(logq(gaussian_filter(f['dm_m384'][:], 1.0)), lo, hi, DMF, 4))
    # inner cube: 12.5 kpc/h voxels vs 23.4 outside -> same mass-density scale needs a log10(volume ratio) shift
    sh_in = 3 * np.log10((9000. / 384) / (2400. / 192))
    wr('idm192.u8.gz', u8(logq(gaussian_filter(f['dm_mI'][:], 1.0)), lo - sh_in, hi - sh_in, DMF, 4))
    wr('shock384.u8.gz', shock_pack(f['max_mach384'][:], f['gas_ed384'][:], S))
    # dissipation is summed over cells: a shock sheet fills a voxel's area, not its volume, so the inner
    # cube needs less than the volume shift; +0.5 dex brings the inner/outer median to agree in the overlap
    wr('ishock192.u8.gz', shock_pack(f['max_machI'][:], f['gas_edI'][:], S, -sh_in + 0.5))

def pack_epoch(snap):
    S = json.load(open(SC)); f = h5py.File(RAW + '/raw_%03d.h5' % snap, 'r'); p = 'ep%03d_' % snap
    print('snap', snap, 'z=%.3f' % f.attrs['z'], 'centre offset from z=0 (ckpc/h):', np.round(f.attrs['cen'] - f.attrs['cen139'], 1),
          'M200=%.2e Msun/h' % (f.attrs.get('group_M200', 0) * 1e10), 'R200=%.0f' % f.attrs.get('group_R200', 0))
    rho = u8(logq(f['gas_m'][:]), S['rho'][0] + LOG8, S['rho'][1] + LOG8)
    L = f['st_L'][:]; sl = u8(logq(L), S['slum'][0] + LOG8, S['slum'][1] + LOG8)
    gr = np.where(L > 0, f['st_Lc'][:] / np.maximum(L, 1e-30), S['scol'][0]); sc = u8(gr, S['scol'][0], S['scol'][1])
    pk = np.stack([rho, sl, sc], axis=-1); wr(p + 'pk192.u8.gz', pk)
    wr(p + 'temp192.u8.gz', u8(logT(f), S['temp'][0], S['temp'][1]))
    wr(p + 'xray192.u8.gz', u8(logq(f['gas_xr'][:]), S['xray'][0], S['xray'][1]))
    wr(p + 'dm192.u8.gz', u8(logq(gaussian_filter(f['dm_m'][:], 1.0)), S['dm'][0] + LOG8, S['dm'][1] + LOG8, DMF, 4))
    wr(p + 'shock192.u8.gz', shock_pack(f['max_mach'][:], f['gas_ed'][:], S, LOG8))
    if 'gas_m384' in f:      # late epochs: density + star light at the z=0 resolution (extract_cubes.py "hires")
        rho = u8(logq(gaussian_filter(f['gas_m384'][:], 1.0)), S['rho'][0], S['rho'][1])
        L = f['st_L384'][:]; sl = u8(logq(L), S['slum'][0], S['slum'][1])
        gr = np.where(L > 0, f['st_Lc384'][:] / np.maximum(L, 1e-30), S['scol'][0]); sc = u8(gr, S['scol'][0], S['scol'][1]); del L, gr
        wr(p + 'pk384.u8.gz', np.stack([rho, sl, sc], axis=-1))
    if 'sp_lr' in f and 'sp_min' in S: wr(p + 'stars.bin.gz', stars_pack(f, S))
    if 'cg_lrho' in f and 'cg_n' in S: wr(p + 'gas.bin.gz', gas_pack(f, S))
    if 'bh_mass' in f or 'lab_pos' in f: meta_json(f, p + 'meta.json')

def galaxies():
    """data/galaxies.json: the subhaloes of the z=0 catalogue inside the cube with M* >= 1e9 Msun, for the page's
    galaxy inspector: p (cube units), ms/mg/mh = log10 of stellar, gas, total subhalo mass [Msun], sfr [Msun/yr],
    rh = stellar half-mass radius [kpc], v = Vmax [km/s], gr = g-r, bh = log10 M_BH, d = distance from the BCG [Mpc]
    (h = 0.681 taken out), i = subhalo index"""
    H = 0.681; cen = np.load(RAW + '/cen139.npy'); box = None; P = []; MT = []; SF = []; RH = []; VM = []; PH = []; BH = []
    for gf in sorted(glob.glob('/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/groups_139/fof_subhalo_tab_139.*.hdf5')):
        with h5py.File(gf, 'r') as g:
            box = float(g['Header'].attrs['BoxSize'])
            if 'Subhalo' not in g or 'SubhaloPos' not in g['Subhalo']: continue
            S = g['Subhalo']; P.append(S['SubhaloPos'][:]); MT.append(S['SubhaloMassType'][:]); SF.append(S['SubhaloSFR'][:])
            RH.append(S['SubhaloHalfmassRadType'][:, 4]); VM.append(S['SubhaloVmax'][:]); PH.append(S['SubhaloStellarPhotometrics'][:]); BH.append(S['SubhaloBHMass'][:])
    P = np.concatenate(P).astype(np.float64); MT = np.concatenate(MT); SF = np.concatenate(SF); RH = np.concatenate(RH); VM = np.concatenate(VM); PH = np.concatenate(PH); BH = np.concatenate(BH)
    d = P - cen; d -= np.round(d / box) * box
    ms = MT[:, 4] * 1e10 / H; keep = np.where((np.max(np.abs(d), axis=1) < HALF) & (ms >= 1e9))[0]
    keep = keep[np.argsort(-ms[keep])]
    out = []
    for i in keep:
        out.append({'i': int(i), 'p': [round(float(v), 4) for v in d[i] / HALF], 'ms': round(float(np.log10(ms[i])), 2),
                    'mg': round(float(np.log10(max(MT[i, 0] * 1e10 / H, 1.))), 2), 'mh': round(float(np.log10(MT[i].sum() * 1e10 / H)), 2),
                    'sfr': round(float(SF[i]), 3), 'rh': round(float(RH[i] / H), 1), 'v': round(float(VM[i]), 0),
                    'gr': round(float(PH[i, 4] - PH[i, 5]), 2), 'bh': round(float(np.log10(max(BH[i] * 1e10 / H, 1.))), 2),
                    'd': round(float(np.linalg.norm(d[i]) / H / 1000.), 3)})
    json.dump(out, open(os.path.join(DATA, 'galaxies.json'), 'w'), separators=(',', ':'))
    print('wrote galaxies.json: %d subhaloes with M* >= 1e9 Msun inside the cube (%.0f kB); M* pct 50/90/max = %s; SFR > 0: %d' % (
        len(out), os.path.getsize(os.path.join(DATA, 'galaxies.json')) / 1e3, [o['ms'] for o in [out[len(out) // 2], out[len(out) // 10], out[0]]], sum(o['sfr'] > 0 for o in out)))

VMAX = 1500.                 # velocity byte: 128 + v / VMAX * 127 (km/s, relative to the hot gas within 500 kpc/h)
SIG = (1.0, 3.0)             # velocity dispersion byte: log10 sigma [km/s] over 10 .. 1000
BV0, BVMAX = 0.05, 10.       # B vector byte: 128 + 127 asinh(B / BV0) / asinh(BVMAX / BV0)  (signed, 0.01 uG still resolved)
MET = (-2.0, 0.5)            # metallicity byte: log10 Z/Zsun
def raw_files():
    return [(int(os.path.basename(q)[4:7]), q) for q in sorted(glob.glob(RAW + '/raw_*.h5'))]

def epname(snap, base):
    return base if snap == 139 else 'ep%03d_' % snap + base

def pack_kin():
    """vel192: RGBA8 = hot-gas velocity vector relative to the cluster's systemic velocity (the mass-weighted mean of
    the hot gas within 500 kpc/h of the centre) and log10 of the in-voxel velocity dispersion (extract_cubes.py kin)"""
    S = json.load(open(SC)); S['vel'] = [VMAX]; S['sig'] = list(SIG); json.dump(S, open(SC, 'w'), indent=1)
    n = 192; c = (np.arange(n) + 0.5) / n * 2 * HALF - HALF; X, Y, Z = np.meshgrid(c, c, c, indexing='ij'); r = np.sqrt(X * X + Y * Y + Z * Z)
    for snap, path in raw_files():
        with h5py.File(path, 'r') as f:
            if 'gas_mhot' not in f: print('snap', snap, 'has no kin part'); continue
            mh = f['gas_mhot'][:]; ok = mh > 0; core = ok & (r < 500.)
            mv = [f['gas_mv' + k][:] for k in 'xyz']; v2 = f['gas_mv2'][:]
        vs = np.array([q[core].sum() / mh[core].sum() for q in mv])
        if snap != 139:   # the epochs at 96^3 (2x2x2 sums, ~2 MB each); z = 0 keeps 192^3 for the XRISM view (19 MB)
            rb = lambda q: q.reshape(96, 2, 96, 2, 96, 2).sum(axis=(1, 3, 5))
            mh = rb(mh); mv = [rb(q) for q in mv]; v2 = rb(v2); ok = mh > 0; n = 96
        else: n = 192
        v = np.stack([np.where(ok, q / np.maximum(mh, 1e-30), 0.) - vs[i] for i, q in enumerate(mv)], -1); v[~ok] = 0.
        s2 = np.where(ok, np.maximum(v2 / np.maximum(mh, 1e-30) - ((v + vs) ** 2).sum(-1), 0.), 0.)
        u = np.zeros((n, n, n, 4), np.uint8)
        # 24 km/s and 0.06 dex steps: the noisy bytes gzip 2-3x better than the full 8 bits
        u[..., :3] = (np.clip(np.round(128 + v / VMAX * 127), 1, 255).astype(np.uint8) // 2) * 2; u[~ok, :3] = 128
        u[..., 3] = u8(0.5 * logq(s2), SIG[0], SIG[1], 0, 8) * ok
        print('snap %d: systemic v = %s km/s, |v| pct 50/90/99 = %s, sigma pct 50/90/99 = %s km/s' % (snap, np.round(vs), np.round(np.percentile(np.linalg.norm(v[ok], axis=1), [50, 90, 99])), np.round(np.sqrt(np.percentile(s2[ok], [50, 90, 99])))))
        wr(epname(snap, 'vel%d.u8.gz' % n), u)

def pack_met():
    S = json.load(open(SC)); S['met'] = list(MET); json.dump(S, open(SC, 'w'), indent=1)
    for snap, path in raw_files():
        with h5py.File(path, 'r') as f:
            if 'gas_mZ' not in f: print('snap', snap, 'has no met part'); continue
            mm = f['gas_mZ_m'][:]; ok = mm > 0; z = np.where(ok, logq(f['gas_mZ'][:] / np.maximum(mm, 1e-30)), -30.)
        print('snap %d: Z/Zsun pct 10/50/90/99 = %s' % (snap, np.round(10 ** np.percentile(z[ok], [10, 50, 90, 99]), 3)))
        wr(epname(snap, 'met192.u8.gz'), u8(z, MET[0], MET[1]))

def pack_bvec():
    S = json.load(open(SC)); S['bvec'] = [BV0, BVMAX]; json.dump(S, open(SC, 'w'), indent=1)
    k = 127. / np.arcsinh(BVMAX / BV0)
    for snap, path in raw_files():
        with h5py.File(path, 'r') as f:
            if 'gas_mBx' not in f: print('snap', snap, 'has no bvec part'); continue
            mm = f['gas_mBv_m'][:]; ok = mm > 0
            B = np.stack([np.where(ok, f['gas_mB' + c][:] / np.maximum(mm, 1e-30), 0.) for c in 'xyz'], -1)
        u = (np.clip(np.round(128 + k * np.arcsinh(B / BV0)), 1, 255).astype(np.uint8) // 2) * 2; u[~ok] = 128
        print('snap %d: |B| pct 50/99 = %s uG, |B_z|/|B| median %.2f' % (snap, np.round(np.percentile(np.linalg.norm(B[ok], axis=1), [50, 99]), 3), np.median(np.abs(B[ok, 2]) / np.maximum(np.linalg.norm(B[ok], axis=1), 1e-9))))
        wr(epname(snap, 'bvec192.u8.gz'), u)

def pack_outer():
    """outer192: RG8 = log total matter mass and log gas mass per 312 kpc/h voxel over +-30 Mpc/h (z=0 only),
    1-voxel smoothed, stretched between the 40th and 99.95th percentiles of the populated voxels"""
    S = json.load(open(SC))
    with h5py.File(RAW + '/raw_139.h5', 'r') as f: am = f['all_mO'][:]; gm = f['gas_mO'][:]
    out = []
    for key, q in (('outer_m', am), ('outer_g', gm)):
        lg = logq(gaussian_filter(q, 1.0)); lo, hi = np.percentile(lg[q > 0], [40, 99.95]); S[key] = [float(lo), float(hi), 1.0]
        print('  %s: lo=%.3f hi=%.3f (populated %d)' % (key, lo, hi, (q > 0).sum())); out.append(u8(lg, lo, hi, 0, 2))
    json.dump(S, open(SC, 'w'), indent=1); wr('outer192.u8.gz', np.stack(out, -1))

def lookback(z, h=0.681, om=0.306):
    """lookback time [Gyr] in flat LCDM"""
    zz = np.linspace(0, z, 2001); E = np.sqrt(om * (1 + zz) ** 3 + 1 - om)
    return float(np.trapz(1. / ((1 + zz) * E), zz) * 977.8 / (100 * h))

def pack_movie():
    """data/mv/mv_SNAP.u8.gz: 128^3 RGB8 (density, log T, star light) of every snapshot on the z=0 scales (the voxel
    is 3x the 384 one: +3 log10 3 on the mass-like fields), T and light kept at 64 levels; data/movie.json: per
    snapshot z, lookback time, M200, R200, centre offset and the tracked label positions (cube units)"""
    S = json.load(open(SC)); sh = 3 * np.log10(3.); os.makedirs(os.path.join(DATA, 'mv'), exist_ok=True); meta = []
    for path in sorted(glob.glob(RAW + '/movie/mv_*.h5')):
        with h5py.File(path, 'r') as f:
            snap = int(f.attrs['snap']); z = max(float(f.attrs['z']), 0.); m = f['gas_m'][:]; lt = np.where(m > 0, f['gas_mlt'][:] / np.maximum(m, 1e-30), 0.)
            L = f['st_L'][:] if 'st_L' in f else np.zeros_like(m)
            # snapshot 139 is centred by cen139.npy, without a group lookup: its M200 / R200 are the page's z = 0 values
            e = {'s': snap, 'z': round(max(z, 0.), 4), 't': round(lookback(max(z, 0.)), 3), 'm': float(f.attrs.get('group_M200', 1.1e5 if snap == 139 else 0)) * 1e10,
                 'r': float(f.attrs.get('group_R200', 1476.35 if snap == 139 else 0)),
                 'off': [round(float(v), 4) for v in (f.attrs['cen'] - f.attrs['cen139']) / HALF]}
            if 'lab_pos' in f:
                names = [q.decode() for q in f.attrs['lab_names']]; lp = f['lab_pos'][:] / HALF; ln = f['lab_n'][:]
                e['lab'] = {nm: [round(float(v), 4) for v in lp[i]] for i, nm in enumerate(names) if np.isfinite(lp[i]).all() and ln[i] >= 20}
        rho = u8(logq(gaussian_filter(m, 0.7)), S['rho'][0] + sh, S['rho'][1] + sh)
        t = u8(lt, S['temp'][0], S['temp'][1], 0, 4); sl = u8(logq(L), S['slum'][0] + sh, S['slum'][1] + sh, 0, 4)
        wr('mv/mv_%03d.u8.gz' % snap, np.stack([rho, t, sl], -1)); meta.append(e)
    meta.sort(key=lambda e: e['s']); movie_flow(meta)
    json.dump(meta, open(os.path.join(DATA, 'movie.json'), 'w'), separators=(',', ':'))
    print('wrote movie.json: %d snapshots, z %.2f .. %.2f, total %.0f MB' % (len(meta), meta[0]['z'], meta[-1]['z'], sum(os.path.getsize(os.path.join(DATA, 'mv', 'mv_%03d.u8.gz' % e['s'])) for e in meta) / 1e6))

def sky():
    """data/sky.json: the observer's frame of this realisation (coma_300/data/zoom_rotations.npz, step 11: the
    observer at the parent box centre, n = Earth -> Coma unit vector in the sim axes, North = +x, the distance from
    coma_300/data/parent_radec.npz) and the real 2M++ galaxies (Lavaux & Hudson 2011, VizieR cone of 6 deg around
    Coma) within 4000 < v_cmb < 10000 km/s as tangent-plane offsets [deg] from NGC 4889 (east +, north +) with Ks"""
    C3 = os.path.expanduser('~/coma_300/data'); STEP = 11
    z = np.load(C3 + '/zoom_rotations.npz'); rd = np.load(C3 + '/parent_radec.npz')
    i = int(np.where(z['step'] == STEP)[0][0]); n = z['n_los'][i]; up = z['up'][i]; east = np.cross(up, n)
    RA0, DEC0 = 195.0338, 27.9770   # NGC 4889
    gal = []
    for line in open(C3 + '/coma_obs/twompp_coma_6deg.tsv'):
        if line.startswith('#') or not line.strip(): continue
        f = line.rstrip('\n').split('\t')
        try: ra, dec, ks, vc = float(f[1]), float(f[2]), float(f[4]), float(f[6])
        except (ValueError, IndexError): continue
        if not 4000 < vc < 10000: continue
        gal.append([round((ra - RA0) * np.cos(np.radians(DEC0)), 4), round(dec - DEC0, 4), round(ks, 2), int(vc)])
    out = {'step': STEP, 'n': [round(float(v), 6) for v in n], 'up': [round(float(v), 6) for v in up], 'east': [round(float(v), 6) for v in east],
           'D': float(rd['dist'][i]), 'ra': float(rd['ra'][i]), 'dec': float(rd['dec'][i]), 'ra0': RA0, 'dec0': DEC0,
           'src': '2M++ (Lavaux & Hudson 2011) within 6 deg of Coma, 4000 < v_cmb < 10000 km/s; frame: coma_300 zoom_rotations step %d' % STEP, 'gal': gal}
    json.dump(out, open(os.path.join(DATA, 'sky.json'), 'w'), separators=(',', ':'))
    print('wrote sky.json: %d galaxies, halo at RA %.2f Dec %.2f, D = %.0f ckpc/h (%.2f deg per Mpc/h)' % (len(gal), out['ra'], out['dec'], out['D'], np.degrees(np.arctan(1000 / out['D']))))

TEMP2 = (5.5, 8.0318)        # temperature byte: log10 T over 3e5 .. 1.1e8 K (the shipped z=0 scale started at 1.3e7 K, which
                             # left 99% of the z = 2 voxels at byte 0 = "cold": the young cluster was tinted cyan and veiled)
def pack_temp():
    """temp192 for z = 0 and every epoch on the TEMP2 scale, all from the raw mass-weighted log T (the z=0 cube from viper
    is replaced: same pipeline at every epoch). Star-forming cells are 1e4 K -> byte 0 = cold, as the shader expects."""
    S = json.load(open(SC)); S['temp_viper'] = S.get('temp_viper', S['temp']); S['temp'] = list(TEMP2) + [1.0]; json.dump(S, open(SC, 'w'), indent=1)
    for snap, path in raw_files():
        with h5py.File(path, 'r') as f: lt = logT(f); m = f['gas_m'][:]
        ok = m > 0; print('snap %d: log T pct 1/10/50/90 = %s, below 1e7 K %.0f%%' % (snap, np.round(np.percentile(lt[ok], [1, 10, 50, 90]), 2), 100 * (lt[ok] < 7).mean()))
        wr(epname(snap, 'temp192.u8.gz'), u8(lt, TEMP2[0], TEMP2[1]))

DMAX = 0.06                  # displacement byte: 128 + d / DMAX * 127, d in cube-width units (9000 ckpc/h) per snapshot interval
def movie_flow(meta):
    """the advected time interpolation: for every snapshot the comoving displacement of the gas to the NEXT snapshot,
    (v - v_core) dt / a on the 48^3 grid (extract_movie.py SNAP more: gas_mV, gas_mv*), as data/mv/dv_SNAP.u8.gz
    (RGB8, DMAX); and data/movie_gal.json: per snapshot the subhaloes [id, x, y, z, log M*, g-r, SF, log M_BH] so
    the page can match neighbours by the most-bound particle id and slide the galaxies along the slider"""
    gal = []; nflow = 0
    for i, e in enumerate(meta):
        path = RAW + '/movie/mv_%03d.h5' % e['s']
        with h5py.File(path, 'r') as f:
            if 'gal' in f:
                g = f['gal'][:]; gal.append([[int(r[0]), round(float(r[1]), 4), round(float(r[2]), 4), round(float(r[3]), 4), round(float(r[4]), 2), round(float(r[5]), 2), int(r[6]), round(float(r[7]), 2)] for r in g])
            else: gal.append([])
            if 'gas_mV' not in f or i == len(meta) - 1: continue
            mV = f['gas_mV'][:]; ok = mV > 0; v = np.stack([np.where(ok, f['gas_mv' + c][:] / np.maximum(mV, 1e-30), 0.) for c in 'xyz'], -1)
        n = mV.shape[0]; c = (np.arange(n) + 0.5) / n * 2 - 1; X_, Y_, Z_ = np.meshgrid(c, c, c, indexing='ij'); core = ok & (np.sqrt(X_ ** 2 + Y_ ** 2 + Z_ ** 2) < 300. / HALF)
        vc = (v[core] * mV[core, None]).sum(0) / mV[core].sum()          # the frame (progenitor) velocity: mass-weighted core mean
        zm = 0.5 * (e['z'] + meta[i + 1]['z']); dt = e['t'] - meta[i + 1]['t']   # Gyr to the next snapshot
        # km/s * Gyr = 1.0227 kpc physical; comoving kpc/h = physical * (1 + z) * h; cube width = 2 HALF ckpc/h
        d = (v - vc) * dt * 1.0227 * (1 + zm) * 0.681 / (2 * HALF); d[~ok] = 0.
        u = np.clip(np.round(128 + d / DMAX * 127), 1, 255).astype(np.uint8); u[~ok] = 128
        wr('mv/dv_%03d.u8.gz' % e['s'], u); e['flow'] = 1; nflow += 1
        if i % 20 == 0: print('   flow %d: |d| pct 50/99 = %s cube widths (max %.3f), dt %.3f Gyr' % (e['s'], np.round(np.percentile(np.linalg.norm(d[ok], axis=1), [50, 99]), 4), np.abs(d).max(), dt))
    json.dump(gal, open(os.path.join(DATA, 'movie_gal.json'), 'w'), separators=(',', ':'))
    print('wrote movie_gal.json (%d snapshots, %d galaxies in all, %.0f kB) and %d flow fields' % (len(gal), sum(len(g) for g in gal), os.path.getsize(os.path.join(DATA, 'movie_gal.json')) / 1e3, nflow))

def pack_mag():
    """magnetic field: mass-weighted |B| [uG] of the non-star-forming gas per 192^3 voxel (extract_cubes.py mag),
    physical uG at every epoch on one log10 scale, 0.001 .. 10 uG (z=0 profile: 3.3 uG within 100 kpc/h, 1.2 at
    300-600, 0.34 at 1-1.5 Mpc/h, 0.11 at 1.5-2.5, 0.017 beyond; peak 11 uG); empty voxels are byte 0"""
    S = json.load(open(SC))
    def lB(f):
        mm = f['gas_mB_m'][:]; return np.where(mm > 0, logq(f['gas_mB'][:] / np.maximum(mm, 1e-30)), -30.), mm > 0
    # 0.001 .. 10 uG: at z = 2 only 6% of the voxels were above a 0.01 uG floor (the young field is weaker), which hid its structure
    S['mag'] = [-3.0, 1.0, 1.0]; json.dump(S, open(SC, 'w'), indent=1); print('wrote', SC, 'mag', S['mag'])
    for path in sorted(glob.glob(RAW + '/raw_*.h5')):
        snap = int(os.path.basename(path)[4:7])
        with h5py.File(path, 'r') as f:
            if 'gas_mB' not in f: print('snap', snap, 'has no gas_mB (run extract_cubes.py %d mag)' % snap); continue
            x, ok = lB(f)
        print('snap %d: %d voxels with ICM gas, |B| [uG] pct 5/50/95/99.9: %s' % (snap, ok.sum(), np.round(10 ** np.percentile(x[ok], [5, 50, 95, 99.9]), 3)))
        wr(('mag192.u8.gz' if snap == 139 else 'ep%03d_mag192.u8.gz' % snap), u8(x, S['mag'][0], S['mag'][1]))

if __name__ == '__main__':
    a = sys.argv[1]
    if a == 'check': sys.exit(0 if check() else 1)
    elif a == 'calib': calib()
    elif a == 'calib_sp': calib_sp()
    elif a == 'calib_cg': calib_cg()
    elif a == 'meta139': meta_json(h5py.File(RAW + '/raw_139.h5', 'r'), 'ep139_meta.json')   # check against D_BH / LABELS, not shipped
    elif a == '139': pack139()
    elif a == 'mag': pack_mag()
    elif a == 'galaxies': galaxies()
    elif a == 'sky': sky()
    elif a == 'temp': pack_temp()
    elif a == 'flow':   # the flow fields + galaxy tracks only (movie.json exists; its `flow` flags are updated)
        meta = json.load(open(os.path.join(DATA, 'movie.json'))); movie_flow(meta); json.dump(meta, open(os.path.join(DATA, 'movie.json'), 'w'), separators=(',', ':'))
    elif a == 'kin': pack_kin()
    elif a == 'met': pack_met()
    elif a == 'bvec': pack_bvec()
    elif a == 'outer': pack_outer()
    elif a == 'movie': pack_movie()
    else: pack_epoch(int(a))
