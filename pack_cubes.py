#!/usr/bin/env python
# Turn the raw cubes of extract_cubes.py into the u8.gz files the page loads.
#   python pack_cubes.py check   -> centre / axis-order check of raw_139.h5 against the shipped rho384
#   python pack_cubes.py calib   -> the check, then fits the u8 scalings of the existing z=0 cubes
#                                   (data/*.u8.gz) against raw_139.h5 and writes cube_scales.json
#   python pack_cubes.py calib_sp -> fits the galaxy-sprite encoding (sp_lum, sp_col, sp_min) of the shipped
#                                   data/stars.bin.gz against the star particles of raw_139.h5 (sp_* datasets)
#   python pack_cubes.py calib_cg -> the cold-gas sprite encoding of data/gas.bin.gz against raw_139.h5 (cg_* datasets)
#   python pack_cubes.py 139     -> data/dm384.u8.gz idm192.u8.gz shock384.u8.gz ishock192.u8.gz
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
    open(os.path.join(DATA, name), 'wb').write(gzip.compress(np.ascontiguousarray(u8).tobytes(), 9))
    print('  wrote', name, '%.1f MB' % (os.path.getsize(os.path.join(DATA, name)) / 1e6), flush=True)

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

def pack_mag():
    """magnetic field: mass-weighted |B| [uG] of the non-star-forming gas per 192^3 voxel (extract_cubes.py mag),
    physical uG at every epoch on one log10 scale, 0.01 .. 10 uG (z=0 profile: 3.3 uG within 100 kpc/h, 1.2 at
    300-600, 0.34 at 1-1.5 Mpc/h, 0.11 at 1.5-2.5, 0.017 beyond; peak 11 uG); empty voxels are byte 0"""
    S = json.load(open(SC))
    def lB(f):
        mm = f['gas_mB_m'][:]; return np.where(mm > 0, logq(f['gas_mB'][:] / np.maximum(mm, 1e-30)), -30.), mm > 0
    if 'mag' not in S:
        S['mag'] = [-2.0, 1.0, 1.0]; json.dump(S, open(SC, 'w'), indent=1); print('wrote', SC, 'mag', S['mag'])
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
    else: pack_epoch(int(a))
