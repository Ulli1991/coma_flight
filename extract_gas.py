import numpy as np, h5py, glob
BASE='/raven/ptmp/uli/sims/borg_coma_zoom_TNG100_first_try/step_011/output/snapdir_139'
with h5py.File('/ptmp/uli/movie011/grids/freeze_139.h5','r') as f:
    cen=f.attrs['cen'][:].astype(np.float64); HALF=float(f.attrs['half_kpch'])
KB=1.3807e-16; MP=1.6726e-24; GAMMA=5./3.
P=[];R=[];T=[];S=[]
for fn in sorted(glob.glob(BASE+'/snap_139.*.hdf5'),key=lambda p:int(p.split('.')[-2])):
    with h5py.File(fn,'r') as f:
        box=f['Header'].attrs['BoxSize']
        c=f['PartType0/Coordinates'][:].astype(np.float64)
        d=c-cen; d-=np.round(d/box)*box
        m=np.max(np.abs(d),axis=1)<HALF
        if not m.any(): continue
        g=f['PartType0']
        rho=g['Density'][:][m]; u=g['InternalEnergy'][:][m]
        xe=g['ElectronAbundance'][:][m]; XH=g['GFM_Metals'][:,0][m]
        sfr=g['StarFormationRate'][:][m]
        mu=4./(1.+3.*XH+4.*XH*xe)
        tt=(GAMMA-1.)*u*1e10*mu*MP/KB
        # ISM / cold dense / star-forming gas only
        sel=(sfr>0)|(tt<10**4.9)
        if not sel.any(): continue
        P.append(d[m][sel].astype(np.float32)); R.append(np.log10(rho[sel]).astype(np.float32))
        T.append(np.log10(np.clip(tt[sel],1e3,None)).astype(np.float32)); S.append((sfr[sel]>0))
pos=np.concatenate(P); lr=np.concatenate(R); lt=np.concatenate(T); sf=np.concatenate(S)
print(len(pos),'cold/SF cells',flush=True)
idx=np.argsort(lr)[::-1][:250000]
np.savez_compressed('/ptmp/uli/movie011/game_gas.npz',
    pos=pos[idx],lrho=lr[idx],ltemp=lt[idx],sf=sf[idx],half=HALF)
print('gas extract done',flush=True)
