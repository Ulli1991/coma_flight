import numpy as np, h5py, glob, gzip
from scipy.spatial import cKDTree
BASE='/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/snapdir_139'
with h5py.File('/ptmp/uli/movie011/grids/freeze_139.h5','r') as f:
    cen=f.attrs['cen'][:].astype(np.float64)
HALF=1200.0; N=192
files=sorted(glob.glob(BASE+'/snap_139.*.hdf5'),key=lambda p:int(p.split('.')[-2]))
P=[];R=[]
for fn in files:
    with h5py.File(fn,'r') as f:
        box=f['Header'].attrs['BoxSize']
        c=f['PartType0/Coordinates'][:].astype(np.float64)
        d=c-cen; d-=np.round(d/box)*box
        m=np.max(np.abs(d),axis=1)<HALF+60.
        if m.any():
            P.append(d[m].astype(np.float32))
            R.append(f['PartType0/Density'][:][m].astype(np.float32))
pos=np.concatenate(P); rho=np.log10(np.clip(np.concatenate(R),1e-12,None))
print(len(pos),'cells',flush=True)
tree=cKDTree(pos)
ax=np.linspace(-HALF,HALF,N,dtype=np.float32)
g=np.empty((N,)*3,np.float32)
for ix in range(N):
    q=np.stack(np.meshgrid(ax[ix:ix+1],ax,ax,indexing='ij'),axis=-1).reshape(-1,3)
    _,ii=tree.query(q,workers=-1)
    g[ix]=rho[ii].reshape(N,N)
lo,hi=np.percentile(g,[40,99.93])
u8=(np.clip((g-lo)/(hi-lo),0,1)*255).astype(np.uint8)
open('/ptmp/uli/movie011/inner_rho.u8.gz','wb').write(gzip.compress(u8.tobytes(),6))
print('inner cube done',flush=True)
