# Dark-matter density cubes for the Coma Virtual Observatory (run on raven/viper, needs the snapshot mount).
# Output: dm384.u8.gz  (384^3, +-HALF kpc/h, log10 density scaled to 0..255)
#         idm192.u8.gz (192^3, +-1200 kpc/h inner core, same scaling recipe)
# Copy both into coma_flight/data/ and push: the page shows the "0 DM" mode when they load.
import numpy as np, h5py, glob, gzip
from scipy.ndimage import gaussian_filter
BASE='/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/snapdir_139'
with h5py.File('/ptmp/uli/movie011/grids/freeze_139.h5','r') as f:
    cen=f.attrs['cen'][:].astype(np.float64); HALF=float(f.attrs['half_kpch'])
N=384; IH=1200.0; NI=192
G=np.zeros((N,)*3,np.float32); GI=np.zeros((NI,)*3,np.float32)
files=sorted(glob.glob(BASE+'/snap_139.*.hdf5'),key=lambda p:int(p.split('.')[-2]))
for fn in files:
    with h5py.File(fn,'r') as f:
        box=f['Header'].attrs['BoxSize']; mt=f['Header'].attrs['MassTable']
        for pt in (1,2):                       # high-res DM, plus low-res boundary DM if present
            k='PartType%d'%pt
            if k not in f or 'Coordinates' not in f[k]: continue
            c=f[k+'/Coordinates'][:].astype(np.float64)
            d=c-cen; d-=np.round(d/box)*box
            m=np.max(np.abs(d),axis=1)<HALF
            if not m.any(): continue
            d=d[m]
            w=f[k+'/Masses'][:][m].astype(np.float64) if 'Masses' in f[k] else np.full(len(d),float(mt[pt]))
            h,_=np.histogramdd(d,bins=N,range=[(-HALF,HALF)]*3,weights=w); G+=h.astype(np.float32)
            mi=np.max(np.abs(d),axis=1)<IH
            if mi.any():
                h,_=np.histogramdd(d[mi],bins=NI,range=[(-IH,IH)]*3,weights=w[mi]); GI+=h.astype(np.float32)
    print(fn,flush=True)
def pack(g,sig,path):
    g=gaussian_filter(g,sig)                  # one-voxel smoothing against particle shot noise
    lg=np.log10(np.clip(g,1e-12,None))
    lo,hi=np.percentile(lg[g>0],[35,99.95])
    u8=(np.clip((lg-lo)/(hi-lo),0,1)*255).astype(np.uint8)
    open(path,'wb').write(gzip.compress(u8.tobytes(),6))
    print(path,'lo/hi',lo,hi,flush=True)
pack(G,1.0,'/ptmp/uli/movie011/dm384.u8.gz')
pack(GI,1.0,'/ptmp/uli/movie011/idm192.u8.gz')
print('dm cubes done',flush=True)
