#!/usr/bin/env python
# Turn the raw cubes of extract_cubes.py into the u8.gz files the page loads.
#   python pack_cubes.py calib   -> fits the u8 scalings of the existing z=0 cubes (data/*.u8.gz)
#                                   against raw_139.h5 and writes cube_scales.json
#   python pack_cubes.py 139     -> data/dm384.u8.gz idm192.u8.gz shock384.u8.gz ishock192.u8.gz
#   python pack_cubes.py 27      -> data/ep027_{pk,xray,temp,dm,shock}192.u8.gz (same scales as z=0)
# All cubes are C-ordered (x,y,z) uint8; the page samples them as uv.zyx.
import sys, os, json, gzip, numpy as np, h5py
from scipy.ndimage import gaussian_filter, maximum_filter
RAW = '/ptmp/uli/coma_cubes'; DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
SC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cube_scales.json')
LOG8 = np.log10(8.)          # a 192^3 voxel holds 8x the mass of a 384^3 voxel
DMF = 28                     # DM shader: smoothstep(0.12,1.0,dd) -> u8 < 31 is invisible

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

def calib():
    f = h5py.File(RAW + '/raw_139.h5', 'r'); S = {}
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
    print('dark matter (384, extract_dm.py recipe: smooth 1, percentiles 35/99.95 of log):')
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

if __name__ == '__main__':
    a = sys.argv[1]
    if a == 'calib': calib()
    elif a == '139': pack139()
    else: pack_epoch(int(a))
