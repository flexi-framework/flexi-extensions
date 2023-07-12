#####################################################################################
# 
# Parametrizes ice shapes read in as series of points in csv format, 
# decomposes the data using principal component analysis (PCA)
# Writes the resulting modes into an HDF5 file.
# 
# Contains lots of plotting routines.
# 
# Ice shapes can be generated from this file by randomly drawing the weighting of
# each mode, calculating a thus weighted sum of modes, and evaluating it at zero 
# to get the new ice outline. 
# 
# Config can be set at beginning of file.
# 
#####################################################################################

import csv 
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import path
from matplotlib import patches
from matplotlib import colors
from matplotlib import cm
from matplotlib.colors import ListedColormap
import sys, os
import copy
import subprocess
import glob
from scipy.interpolate import CubicSpline
import matplotlib._cntr as cntr
from sklearn.decomposition import PCA
import h5py
import commons


#----------------------------------------------------------------------------------------------------
# # Used prms for 8 modes 
min_spray_time = 0.
grid_spacing = 5.E-4
spline_spacing = grid_spacing/2.
xylim = (-0.1,0.3,-0.1,0.1)
maxdist = 0.05
fac_signed_dist_func = 3.
do_smooth_post = False
# do_show = [False, # Orginial Shapes
           # False, # mode contours
           # False, # modes (field) 
           # False, # random samples
           # False, # random sample modes
           # False, # original shape fields
           # False] # colormap
do_show = [True, # Orginial Shapes
           True, # mode contours
           True, # modes (field) 
           True, # random samples
           True, # random sample modes
           True, # original shape fields
           True] # colormap
# savefig_formats = ['jpg','pgf','png']
savefig_formats = ['pgf','png']
n_write_lineplots = 10
do_save_modes = False
do_save_csv = True


#----------------------------------------------------------------------------------------------------

path_afl = "report/airfoil.csv"
x_afl = commons.read_coords(path_afl)
ind = np.argwhere(x_afl[:,0] < 0.25)
# cut profile at 0.25 
x_afl = x_afl[ind[0][0]:ind[-1][0]+1,:]
x_afl = np.flip(x_afl,axis=0)


#----------------------------------------------------------------------------------------------------

filename = "report/metadata.csv"
with open(filename) as csvfile:
    reader = csv.reader(csvfile)
    data = [[e for e in r] for r in reader]
    title = data.pop(0)

names_ice = []
files = glob.glob("report/shapes/*.csv")
# files = ["json/613.csv"]
for f in files: 
    name = f.split("/")[-1].split(".")[0]
    is_long = False
    for r in data: 
        if r[0] == name: 
            is_long = float(r[-1]) > min_spray_time
            break
    if is_long: 
        names_ice.append(name)

#----------------------------------------------------------------------------------------------------


xmin, xmax, ymin, ymax = xylim
nx, ny = int((xmax-xmin)/grid_spacing)+1, int((ymax-ymin)/grid_spacing)+1 
npts = nx*ny

x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)
X,Y = np.meshgrid(x,y)
V = np.moveaxis(np.array([X,Y]),0,-1).reshape((npts,2))

n_data = len(names_ice)
avg = np.zeros((npts,))
al = np.ndarray((n_data,npts))

p_afl = path.Path(x_afl)
c_afl = p_afl.contains_points(V).astype(int)-1.+maxdist

n_xy = int(np.ceil(np.sqrt(len(names_ice))))
if do_show[0]: fig0,ax0=plt.subplots(n_xy,n_xy,sharex=True,sharey=True,figsize=(25,15))
if do_show[1]: fig1,ax1=plt.subplots(n_xy,n_xy,sharex=True,sharey=True,figsize=(25,15))
if do_show[2]: fig2,ax2=plt.subplots(n_xy,n_xy,sharex=True,sharey=True,figsize=(25,15))
if do_show[3]: fig3,ax3=plt.subplots(n_xy,n_xy,sharex=True,sharey=True,figsize=(25,15))
if do_show[4]: fig4,ax4=plt.subplots(n_xy,n_xy,sharex=True,sharey=True,figsize=(25,15))
if do_show[5]: fig5,ax5=plt.subplots(n_xy,n_xy,sharex=True,sharey=True,figsize=(25,15))
if do_show[6]: fig6,ax6=plt.subplots(1,1,figsize=(5,2))

# colormap
mynorm = colors.Normalize(vmin=-1., vmax=1., clip=False)

cmap = plt.get_cmap('viridis')
c1 = cmap(0.15)
c2 = cmap(0.7)
linecol = [0.8*c for c in cmap(0.8)[:3]]
br_rat = np.sqrt(np.sum(c1[:3])/np.sum(c2[:3]))
c1 = [i/br_rat for i in c1]
c2 = [i*br_rat for i in c2]
c3 = cmap(1.)

N = 256
vals = np.ones((N, 4))
for i in range(3): 
    mc = 1. #c3[i]**0.5
    s = np.linspace(0.,1.,128)**1.7
    ti = (1.-s      )*0.8*c1[i] + s      *mc
    bi = (1.-s[::-1])*0.8*c2[i] + s[::-1]*mc
    vals[:,i] = np.concatenate((ti,bi))
newcmp = ListedColormap(vals,name='new')

# alternative colormap
# newcmp = plt.get_cmap('seismic')
# linecol = [0.8, 0., 0.]

#----------------------------------------------------------------------------------------------------

x_write = []
for i_prf, name_ice in enumerate(names_ice): 
# for name_ice in names_ice[0:1]: 
    path_ice = "report/shapes/"+name_ice+".csv"
    x_ice = commons.read_coords(path_ice)
    # add [x, 0.] on each side to make sure we start inside the airfoil
    # x_ice = np.concatenate((np.array([[x_ice[0,0],0.]]),x_ice,np.array([[x_ice[-1,0],0.]])))
    x_ice = np.concatenate((np.array([[0.2,0.03]]),x_ice,np.array([[0.2,0.]])))

    x_out = commons.path_union(x_afl,x_ice)


    # interpolate
    x_out = commons.interpolate(x_out,spac=spline_spacing)
    x_write.append(x_out[x_out[:,0]<0.2,:])

    if do_show[0]:
        ix = i_prf%n_xy
        iy = i_prf//n_xy
        ax0[iy,ix].plot(x_out[:,0],x_out[:,1],color=linecol,lw=2)
        ax0[iy,ix].plot(x_afl[:,0],x_afl[:,1],'k',lw=2)

    p = path.Path(x_out)
    # c = p.contains_points(V,radius=-0.03)
    c = 2.*p.contains_points(V).astype(int)-1.
    d = maxdist*np.ones((ny,nx))


    for v in x_out: 
        il = max(int((v[0]-maxdist-xmin)/grid_spacing)-1,0)
        iu = min(int((v[0]+maxdist-xmin)/grid_spacing)+1,nx)
        jl = max(int((v[1]-maxdist-ymin)/grid_spacing)-1,0)
        ju = min(int((v[1]+maxdist-ymin)/grid_spacing)+1,ny)
        l = np.sqrt((X[jl:ju,il:iu]-v[0])**2+(Y[jl:ju,il:iu]-v[1])**2)
        d[jl:ju,il:iu] = np.minimum(d[jl:ju,il:iu], l)
    c = c*np.reshape(d,(npts,))
    # c = np.maximum(c,c_afl)
    c = np.arctan(fac_signed_dist_func*c/maxdist)
        
    if do_show[5]:
        ix = i_prf%n_xy
        iy = i_prf//n_xy
        tmp = c.reshape((ny,nx))
        ax5[iy,ix].imshow(tmp,extent=xylim,interpolation='bilinear',origin='lower',norm=mynorm,cmap = newcmp)

    al[i_prf,:] = c
    avg += c
avg /= n_data

for i in range(n_data):
    al[i,:] -= avg

pca = PCA()
pca.fit(al)

if do_save_csv: 
    np.savetxt("fig_output/singular_values.csv",pca.singular_values_)

def smooth(t,a,n):
    for j in range(n): 
        t[1:-1,:] = (1.-a)*t[1:-1,:] + a/2.*(t[:-2,:]+t[2:,:])
        t[:,1:-1] = (1.-a)*t[:,1:-1] + a/2.*(t[:,:-2]+t[:,2:])
    return t



def get_cntr(arr):
    c = cntr.Cntr(X, Y, arr)
    nlist = c.trace(0., 0., 0)
    return nlist[:len(nlist)//2][0]

tmp = avg.reshape((ny,nx))
avg = smooth(tmp,0.4,2)
if do_show[2]: 
    ax2[0,0].imshow(avg,extent=xylim,interpolation='bilinear',origin='lower',norm=mynorm,cmap = newcmp)
modes=[]
for i in range(n_data-1):
    ix = (i+1) % n_xy
    iy = (i+1)//n_xy
    tmp = pca.singular_values_[i]*pca.components_[i,:]
    tmp2 = tmp.reshape((ny,nx))
    tmp = smooth(tmp2,0.4,2)
    modes.append(tmp)
    if do_show[1]: 
        segs= get_cntr(avg+tmp)
        ax1[iy,ix].plot(segs[:,0],segs[:,1],color=linecol,lw=2)
        segs= get_cntr(avg-tmp)
        ax1[iy,ix].plot(segs[:,0],segs[:,1],color=linecol,lw=2)
        ax1[iy,ix].plot(x_afl[:,0],x_afl[:,1],'k',lw=2)
    if do_show[2]: 
        ax2[iy,ix].imshow(tmp,extent=xylim,interpolation='bilinear',origin='lower',norm=mynorm,cmap = newcmp)
if do_show[1]: 
    segs = get_cntr(avg)
    ax1[0,0].plot(segs[:,0],segs[:,1],color=linecol,lw=2)
    ax1[0,0].plot(x_afl[:,0],x_afl[:,1],'k',lw=2)

def get_sample(i_sample): 
    # np.random.seed(i_sample+70)
    # # a = np.sqrt(3.)/2.
    a = 0.75
    # phi = np.random.uniform(-a,a,n_data-1)
    # # phi = np.random.normal(0.,0.5,n_data-1)

    tmp = copy.deepcopy(avg)
    # for phii,modei in zip(phi,modes): 
    for i_mode,mode in enumerate(modes): 
             #          seed_id
        seed = int(1e8)*0 + int(1e2)*(i_sample) + i_mode + 50
        np.random.seed(seed)
        phi = np.random.uniform(-a,a,1)
        tmp += phi*mode
    if do_smooth_post: 
        tmp = smooth(tmp,0.4,10)
    return tmp

samples_write = []
if do_show[3]: 
    for i_sample in range(n_data):
        ix = i_sample % n_xy
        iy = i_sample//n_xy

        tmp = get_sample(i_sample)
        segs= get_cntr(tmp)

        samples_write.append(segs)
        ax3[iy,ix].plot(segs[:,0],segs[:,1],color=linecol,lw=2)
        ax3[iy,ix].plot(x_afl[:,0],x_afl[:,1],'k',lw=2)

if do_show[4]: 
    for i_sample in range(n_data):
        ix = i_sample % n_xy
        iy = i_sample//n_xy

        tmp = get_sample(i_sample)

        ax4[iy,ix].imshow(tmp,extent=xylim,interpolation='bilinear',origin='lower',norm=mynorm,cmap = newcmp)

if do_show[6]: 
    d1 = np.linspace(-1.,1.,N)
    d2 = np.linspace(0.,0.1,2)
    dm,_ = np.meshgrid(d1,d2)
    ax6.imshow(dm,extent=(0.,1.,0.,0.1),interpolation='bilinear',origin='lower',norm=mynorm,cmap = newcmp)


if do_show[0]: 
    ax0[0,0].axis('equal')
    ax0[0,0].axis((-0.05,0.15,-0.06,0.08))
    ax0[0,0].set_xticks([])
    ax0[0,0].set_yticks([])
    for fmt in savefig_formats:
        fig0.savefig('fig_output/original_shapes.'+fmt)
if do_show[1]: 
    ax1[0,0].axis('equal')
    ax1[0,0].axis((-0.05,0.15,-0.06,0.08))
    ax1[0,0].set_xticks([])
    ax1[0,0].set_yticks([])
    for fmt in savefig_formats:
        fig1.savefig('fig_output/mode_contours.'+fmt)
if do_show[2]: 
    ax2[0,0].axis('equal')
    ax2[0,0].axis((-0.05,0.15,-0.06,0.08))
    ax2[0,0].set_xticks([])
    ax2[0,0].set_yticks([])
    for fmt in savefig_formats:
        fig2.savefig('fig_output/mode_fields.'+fmt)
if do_show[3]: 
    ax3[0,0].axis('equal')
    ax3[0,0].axis((-0.05,0.15,-0.06,0.08))
    ax3[0,0].set_xticks([])
    ax3[0,0].set_yticks([])
    for fmt in savefig_formats:
        fig3.savefig('fig_output/random_samples.'+fmt)
if do_show[4]: 
    ax4[0,0].axis('equal')
    ax4[0,0].axis((-0.05,0.15,-0.06,0.08))
    ax4[0,0].set_xticks([])
    ax4[0,0].set_yticks([])
    for fmt in savefig_formats:
        fig4.savefig('fig_output/sample_fields.'+fmt)
if do_show[5]: 
    ax5[0,0].axis('equal')
    ax5[0,0].axis((-0.05,0.15,-0.06,0.08))
    ax5[0,0].set_xticks([])
    ax5[0,0].set_yticks([])
    for fmt in savefig_formats:
        fig5.savefig('fig_output/original_shape_fields.'+fmt)
if do_show[6]: 
    ax6.axis('equal')
    ax6.axis((0.,1.,0.,0.1))
    ax6.set_xticks([])
    ax6.set_yticks([])
    for fmt in savefig_formats:
        fig6.savefig('fig_output/colormap.'+fmt)

if do_save_modes: 
    modes = np.array(modes)
    filename='modes.h5'
    with h5py.File(filename, 'w') as h5f:
         h5f.create_dataset('average', data=avg)
         h5f.create_dataset('modes', data=modes)
         h5f.create_dataset('x', data=x)
         h5f.create_dataset('y', data=y)

if do_save_csv: 
    for i in range(n_data): 
        np.savetxt("fig_output/x_afl.csv",x_afl[x_afl[:,0]<0.2,:])
        np.savetxt("fig_output/original_shapes_"+str(i)+".csv",x_write[i][::n_write_lineplots])
        np.savetxt("fig_output/samples_"+str(i)+".csv",samples_write[i][::n_write_lineplots])



plt.show()





# patch = patches.PathPatch(p, facecolor='none', lw=2)
# ax.add_patch(patch)






