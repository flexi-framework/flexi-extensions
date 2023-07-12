#####################################################################################
# 
# Generates hybrid meshes around iced aifroil shapes
# Consiting of structured boundary layer meshes distorted with an in-house algorithm
# and an unstructred farfield mesh. 
# 
# Requires modes file created with get_pca_modes.py script. 
# 
# For use, see arguments at end of script. 
# 
#####################################################################################

import h5py 
import numpy as np
import copy
import sys
from matplotlib import cm
import matplotlib._cntr as cntr
import csv 
from scipy.interpolate import CubicSpline,interp1d,interp2d,UnivariateSpline
import glob 
from matplotlib.collections import LineCollection
import meshio 
import pygmsh
import warnings
import subprocess

import commons


def make_mesh(hopr_path,reference_mode,hoprbasefile,name_str,random_vars,n_avg):

   # --------------------------------------------------------------------------------
   # INPUT
   # --------------------------------------------------------------------------------
   # FIGURES (do not turn on on hawk)
   do_show_fig = False
   # do_show_fig = True
   do_plot_inner_lines = True
   inner_line_str = ":"
   n_steps_show = 200

   ax_limits=(-0.1,0.1,-0.1,0.1)
   # ax_limits=(-0.05,0.05,-0.05,0.05)
   # ax_limits=(-0.1,0.1,-0.1,0.1)
   # ax_limits=(-0.9,0.9,-0.9,0.9)

   # --------------------------------------------------------------------------------
   # output 
   do_save_fig = False

   do_write_mesh = True
   # do_write_mesh = False

   # struct mesh prms
   n_geo = 3

   # OLD - determine via spacing instead
   # n_elems_front = 40
   # n_front = n_elems_front*n_geo+1


   x_max_front = 0.1                 # pos where streamwise mesh stretching starts
   spac_front = 0.0015                   # streamwise mesh spacing at x_max_front
   spac_int = 0.004                   # max streamwise mesh spacing on profile
   spac_back = 0.0025                # streamwise mesh stretching at trailing edge
   stretch_int = 1.02                 # stretch factor for main profile

   spacmax_fac = 3.5                    # spac_max / spac_front
   spacmax = spac_front * spacmax_fac
   scale_fac_front_ref = 1.2 # iced lines can be longer than clean ones. This factor makes front spacing smaller for clean airfoil in reference mode.

   drad0 = 0.00035                     # first layer thickness
   n_elems_rad = 25       
   nrad = n_elems_rad*n_geo + 1
   stretch_fac = 1.05

   # re-sample
   post_drad0 = 0.00025
   post_dradmax = 0.004
   post_stretch = 1.07
   post_n_elems_rad = 18
   post_nrad = post_n_elems_rad * n_geo + 1

   # prms and constraints for movement
   n_steps = 2000
   n_steps_after = 1000

   dpamin = 1.1*drad0

   alpha_max = 60*np.pi/180
   laplace_amplitude = 0.4
   max_ratio_move_fac = 0.6
   max_ratio = 1.5

   interp_spacing = 2E-4       # determines interpoaltion accuracy, not final spacing

   do_remove_inner_stretching = True

   # unstructured mesh 
   r_circ = 10.
   circ_ctr = [1.0, 0., 0.]
   al_wake = 3.*np.pi/180.           # AoA 
   phi_wake =12.*np.pi/180.          # half widening angle: min/max = alpha +- phi
   n_circ = 52                       # Ice
   # n_circ = 48                       # no Ice
   stretch_fac_wake = 1.035 

   # ice shape
   do_ice = True

   # --------------------------------------------------------------------------------
   # --------------------------------------------------------------------------------


   if do_show_fig or do_save_fig:
       plt.ion()
       plt.figure(figsize = (10,10))
       fig = plt.gcf()

   min_ = np.minimum
   max_ = np.maximum

   def constrain_0to1(arr): 
       return min_(1.,max_(0.,arr))

   def length(arr): 
       return np.linalg.norm(arr,axis = -1)

   def Laplace_inner(X,gamma): 
       return (1.-gamma)*X[2:-1,1:-1,:] + gamma*0.25*(X[3:,1:-1,:] + X[1:-2,1:-1,:] + X[2:-1,2:,:] + X[2:-1,:-2,:])

   def angle(a,b): 
       la = length(a)
       lb = length(b)
       c = np.sum(a*b,axis=-1)/(la*lb)
       c = max_(min_(c,1.),-1.)
       return np.arccos(c), la, lb
       # <a,b> = |a|*|b|*cos(alpha)

   def add_dim(arr): 
       return arr.reshape(arr.shape + (1,))


   # FILENAMES
   filename_coords = name_str+"_coords2D.h5"
   basename = "_".join(name_str.split("_")[:-3])
   filename_coords_base = basename+"_coords2D.h5"
   filename_geo = basename+"_unstruct.geo"
   filename_msh = basename+"_unstruct.msh"

   # --------------------------------------------------------------------------------
   # READ CLEAN AIRFOIL COORDS

   path_afl = "report/airfoil.csv"
   x_afl = commons.read_coords(path_afl)

   # SHARP LEADING EDGE
   # extrapolate
   y_last = x_afl[-1,1]+(x_afl[-1,1]-x_afl[-2,1])/(x_afl[-1,0]-x_afl[-2,0])*(1.-x_afl[-1,0])
   x_afl = np.concatenate((x_afl,np.array([[1.,y_last]])))
   # adjust suction side
   dy = x_afl[0,1]-y_last
   h = x_afl.shape[0]//2 # leading edge int (roughly)
   # sine ramping between x=0.3 and x=1.0
   fac = 0.5-0.5*np.cos(np.pi*constrain_0to1((x_afl[h:,0]-0.3)/0.7))
   x_afl[h:,1]+=fac*dy

   # --------------------------------------------------------------------------------
   # GET ICE SHAPE FROM MODES 

   # read in modes
   filename='modes.h5'
   with h5py.File(filename, 'r') as h5f:
           avg = np.array(h5f['average'])
           modes = np.array(h5f['modes'])
           x = np.array(h5f['x'])
           y = np.array(h5f['y'])
   X,Y = np.meshgrid(x,y)

   assert len(random_vars) <= len(modes)

   # linear combination
   tmp = copy.deepcopy(avg)
   for phii,modei in zip(random_vars,modes): 
       tmp += float(phii)*modei

   # isoline / contour is ice outline
   c = cntr.Cntr(X, Y, tmp)
   nlist = c.trace(0., 0., 0)
   segs = nlist[:len(nlist)//2][0]

   # --------------------------------------------------------------------------------
   # JOIN ICE SHAPE AND CLEAN AIRFOIL

   # only use left (front) part of closed outline 
   x_cut = 0.21
   n = int(0.04/interp_spacing)
   if segs[0,0] > x_cut: 
       i_use = np.argwhere(segs[:,0]<x_cut)
       il = i_use[0][0]
       iu = i_use[-1][0] 
       # add first and last points inside airfoil
       # xl = [0.25, segs[il,1]]
       # xu = [0.25, segs[iu,1]]
       xl = np.transpose(np.array([np.linspace(0.25,segs[il,0]+interp_spacing,n), segs[il,1]*np.ones((n,))]))
       xu = np.transpose(np.array([np.linspace(segs[iu,0]+interp_spacing,0.25,n), segs[iu,1]*np.ones((n,))]))
       segs = np.concatenate((xl, segs[il:iu+1,:], xu))
   else: 
       i_discard = np.argwhere(segs[:,0]>x_cut)
       il = i_discard[0][0]-1
       iu = i_discard[-1][0]+1
       # add first and last points inside airfoil
       # xl = [0.25, segs[il,1]]
       # xu = [0.25, segs[iu,1]]
       xl = np.transpose(np.array([np.linspace(segs[il,0]+interp_spacing,0.25,n), segs[il,1]*np.ones((n,))]))
       xu = np.transpose(np.array([np.linspace(0.25,segs[iu,0]+interp_spacing,n), segs[iu,1]*np.ones((n,))]))
       segs = np.concatenate((xu, segs[iu:,:], segs[:il+1,:], xl))
   if segs[0,1] > segs[-1,1]: 
       segs = np.flip(segs,axis=0)

   # inperpolate both clean airfoil and ice outline  
   segs = commons.interpolate(segs,spac=interp_spacing)
   alpha=0.5
   for i in range(10):
       segs[1:-1,:] = (1.-alpha)*segs[1:-1,:] + alpha*0.5*(segs[2:,:]+segs[:-2,:])
   x_afl_n = commons.interpolate(x_afl,spac=interp_spacing)
   # use outermost of the two (ice thickness can't be negative)
   segs = commons.path_union(x_afl_n,segs)

   # l = commons.linlen(segs)
   # spx = UnivariateSpline(l, segs[:,0],k=3,s=1.e-4)
   # spy = UnivariateSpline(l, segs[:,1],k=3,s=1.e-4)


   # --------------------------------------------------------------------------------
   # INTERPOLATE ACCORDING TO SPACING ALONG CLEAN AND ICED AIRFOIL OUTLINE

   # procedure (for each): 
   # - get lentgh of polyline 
   # - cut ar x_max_front to divide into lower back, front, and upper back 
   # - get l positions for each of the three
   #    - spacing for back is stretching
   #    - spacing for front is based on concaveness (see details below)
   # - join lower back, fron and upper back to one
   # - get spline through points and interpolate to new pos

   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   # functions 

   # get linear length of polyline and length at x = x_max_front
   def get_l_lr(arr): 
       l = commons.linlen(arr)
       xsp = CubicSpline(l,arr[:,0]-x_max_front)
       lr = xsp.roots()
       lr = [i for i in lr if 0. < i < l[-1]]
       if len(lr) != 2: 
           print(lr)
           sys.exit()
       return l, lr

   # get y with constant stretching with inputs:
   # overall length, first spacing, last spacing
   def get_distrib(l_overall,l_first,l_last,l_int,stretch,force_n=0,uneven=-1): 
       nf = int(np.log(l_int/l_first)/np.log(stretch))+1
       af = np.array([l_first*stretch**i for i in range(nf)])
       lf = np.sum(af)
       nb = int(np.log(l_int/l_last)/np.log(stretch))+1
       ab = np.array([l_last*stretch**i for i in range(nb)])
       lb = np.sum(ab)
       if force_n: 
           n = force_n - 1 # leading zero is added
       else: 
           n = int(nf+nb+l_overall/l_int)+1
           if n % n_geo: 
               n += n_geo - (n % n_geo)
           if uneven >= 0: # -1 means deactivated
               n += ((uneven+n//n_geo)%2) * n_geo
       ni = n-nf-nb
       dxs = np.concatenate((af,l_int*np.ones((ni,)),np.flip(ab)))
       xs = np.concatenate((np.array([0.]),np.cumsum(dxs)))
       xs = xs*l_overall/xs[-1]
       return xs



   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   # iced airfoil

   # length of polyline  
   l_iced, l_cut_iced = get_l_lr(segs)
   # get spline 
   sp=CubicSpline(l_iced,segs)

   # if not reference mode: force n pts in each segment
   if reference_mode: 
       force_n_upper = 0
       force_n_lower = 0
   else: 
       with h5py.File(filename_coords_base, 'r') as h5f:
           force_n_upper = h5f.attrs['n_ref_upper']
           force_n_front = h5f.attrs['n_ref_front']
           force_n_lower = h5f.attrs['n_ref_lower']

   # spacing in iced front is larger in regions with concave ice shapes. 
   # strategy: base strecthing on size of circle that would still be able to touch this region from outside 

   # get spline for outline angle
   sp_x   = CubicSpline(l_iced,segs[:,0])
   sp_y = CubicSpline(l_iced,segs[:,1])
   sp_xp  = sp_x.derivative(1)
   sp_yp  = sp_y.derivative(1)
   alpha = np.arctan2(sp_xp(l_iced),sp_yp(l_iced))
   for i in range(alpha.shape[0]-1):
       if alpha[i+1]-alpha[i] > np.pi: 
           alpha[i+1:]=alpha[i+1:]-2*np.pi
       elif alpha[i+1]-alpha[i] < -np.pi: 
           alpha[i+1:]=alpha[i+1:]+2*np.pi
   sp_alpha  = CubicSpline(l_iced,alpha)

   # new splines (just the front part)
   nf = int((l_cut_iced[1]-l_cut_iced[0])/interp_spacing)
   lf = np.linspace(l_cut_iced[0], l_cut_iced[1], nf) 
   sn = sp(lf)
   an = sp_alpha(lf)



   # ROADMAP: 
       # - get r/ds in matrix form
       # - find max criterion based on dpa_min
       # - average until small enough 
       # - adapt sn and lf (no need for new splines) 

   warnings.filterwarnings("ignore", category=RuntimeWarning)

   # first index is loop index = row
   eps = 1.e-6
   x,y = np.meshgrid(sn[:,0],sn[:,1],indexing='ij')
   anm,_ = np.meshgrid(an,an,indexing='ij')
   dx = np.transpose(x)-x
   dy = y-np.transpose(y)
   den = 2.*(-dx*np.cos(anm)+dy*np.sin(anm))
   enu = dx**2 + dy**2 
   tmp = np.where(den>eps,enu/den,10.*np.ones_like(anm))
   minvalv = np.amin(tmp,axis=1)
   minvalm, _ = np.meshgrid(minvalv,minvalv,indexing='ij')
   jv = np.argmin(tmp,axis=1)
   iv = np.arange(nf) 
   iminv = min_(iv,jv)
   imaxv = max_(iv,jv)
   lv = (imaxv-iminv)*interp_spacing
   iminm, imaxm = np.meshgrid(iminv,imaxv,indexing='ij')
   imaxm = np.transpose(imaxm)
   _,im = np.meshgrid(iv,iv,indexing='ij')
   mask = np.logical_and( iminm <= im , im <= imaxm)
   r = np.where(mask,minvalm,10.)
   r = np.amin(r,axis=0)

   exc_max = 0.25
   n_average = 500
   n_transition = 20

   npa = (imaxv-iminv)*interp_spacing/spacmax
   dpa_exc = (npa*dpamin)/(2*r) #should be < 1 
   do_average = dpa_exc>exc_max
   if np.sum(do_average): 
       do_average, _ = np.meshgrid(do_average,do_average,indexing='ij')
       do_average = np.logical_and( mask,do_average>0.)
       do_average = np.amax(do_average,axis=0)
       do_average = add_dim(do_average)[1:-1,:]
       # sn_orig = copy.deepcopy(sn)
       for i in range(n_average): 
           if i+n_transition>n_average: 
               do_average[1:-1] = (do_average[2:]+do_average[1:-1]+do_average[:-2])>0.
           sn[1:-1,:] = (1-do_average)*sn[1:-1,:] + do_average*0.5*(sn[2:,:]+sn[:-2,:])
       lf = commons.linlen(sn)
       lf = l_cut_iced[0] + (l_cut_iced[1]-l_cut_iced[0])*lf/lf[-1]

   # fig = plt.gcf()
   # fig.clf()
   # ax = plt.gca()
   # ax.axis('equal')
   # ax.set_xticks([])
   # ax.set_yticks([])
   # ax.axis((-0.1,0.1,-0.1,0.1))
   # plt.plot(sn[:,0],sn[:,1],'r')
   # plt.plot(sn_orig[:,0],sn_orig[:,1],'-+k')
   # # plt.scatter(sn_orig[:,0],sn_orig[:,1],c=cm.hot((minvalv-np.amin(minvalv)-0.005)/-0.005), edgecolor='none')
   # fig.canvas.draw_idle()  # required to work on mpl < 1.5
   # fig.canvas.flush_events()
   # input()
   # sys.exit()



   tmp2 = np.where(den<-eps,-enu/den,10.*np.ones_like(anm))
   minvalv = np.amin(tmp2,axis=1)
   minvalm, _ = np.meshgrid(minvalv,minvalv,indexing='ij')
   jv = np.argmin(tmp2,axis=1)
   iv = np.arange(nf) 
   iminv = min_(iv,jv)
   imaxv = max_(iv,jv)
   iminv, imaxv = np.meshgrid(iminv,imaxv,indexing='ij')
   imaxv = np.transpose(imaxv)
   _,im = np.meshgrid(iv,iv,indexing='ij')
   mask = np.logical_and( iminv <= im , im <= imaxv)
   r_inner = np.where(mask,minvalm,10.)
   r_inner = np.amin(r,axis=0)

   warnings.filterwarnings("default", category=RuntimeWarning)

   # l = commons.linlen(segs)
   # l1, l2 = np.meshgrid(l,indexing='ij')
   # dl = np.abs(l1-l2)

   # for each point, calculate which size circle is still able to touch it
   # eps = 1.e-6
   # r = 10.*np.ones_like(lf) # init circle radius with large value
   # r_inner = 10.*np.ones_like(lf) # init circle radius with large value
   # # runtime warning disabled for divide by zero (enu=0.), which is discarded in np.where
   # warnings.filterwarnings("ignore", category=RuntimeWarning)
   # for i in range(nf): 
       # # distance to all other points 
       # dx = sn[:,0]-sn[i,0]
       # dy = sn[:,1]-sn[i,1]
       # # denominator = 2.* < dx , n > with dx vector to other points and n unit surf normal vector
       # den = 2.*(-dx*np.cos(an[i])+dy*np.sin(an[i]))
       # # enumerator = |dx|**2
       # enu = dx**2 + dy**2
       # # circle radius, if it is positive, else large value
       # tmp = np.where(den>eps,enu/den,10.*np.ones_like(lf))
       # # diminish r between two touch points of circle
       # minval = np.min(tmp)
       # j = np.argmin(tmp)
       # imin, imax = min(i,j), max(i,j)+1
       # r[imin:imax] = min_(r[imin:imax],minval)

       # # circle radius, if it is positive, else large value
       # tmp = np.where(den<-eps,-enu/den,10.*np.ones_like(lf))
       # # diminish r locally
       # minval = np.min(tmp)
       # j = np.argmin(tmp)
       # r_inner[i] = min_(r_inner[i],minval)
       # r_inner[j] = min_(r_inner[j],minval)
       
   # print(np.sum(r))
   # print(r)
   warnings.filterwarnings("default", category=RuntimeWarning)
   # sys.exit()



   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   # base spacing on r: weighted average (large exponent increases importance of small r)
   expo = 10.
   spac = r**(1./expo)
   alpha = 1.
   for i in range(30): 
       spac[1:-1] = (1.-alpha)*spac[1:-1] + alpha*0.5*(spac[2:]+spac[:-2])
   curve = spac**-expo

   # adjust to min and max spacing, cap
   max_abs_curv = 250.
   min_abs_curv = 0.
   spac = spac_front + (spacmax-spac_front) * (max_(min_(curve,max_abs_curv),min_abs_curv+eps)-min_abs_curv)/(max_abs_curv-min_abs_curv)

   # enlarge concave regions with large spacing
   for i in range(5): 
       spac[1:-1] = max_(max_(spac[1:-1], spac[2:]), spac[:-2])

   # smaller spacing in highly curved regions
   spacmin = spac_front/1.5
   max_abs_curv = 600.
   min_abs_curv = 300.
   spac_inner = spacmax + (spacmin-spacmax) * constrain_0to1((1./r_inner-min_abs_curv)/(max_abs_curv-min_abs_curv))
   spac = min_(spac,spac_inner)

   # get actual coordinates from theoretical spacing values at equistant coordinates
   dens = 1./spac
   intdens = np.cumsum((lf[1:]-lf[:-1])*0.5*(dens[1:]+dens[:-1]))
   xi = np.concatenate((np.array([0.]),intdens))
   xi_l = 0.
   xi_u = xi[-1]

   # number of pts in front part of airfoil
   if reference_mode: 
       n_front_float = (xi_u - xi_l)*scale_fac_front_ref
       n_front = int(np.ceil(n_front_float/n_geo)*n_geo)+1
   else: 
       n_front = force_n_front

   sp_segs_of_xi = CubicSpline(xi,sn)
   # unit spacing in reference space 
   xi_eval = np.linspace(xi_l,xi_u,n_front)
   # transform to physical space
   xsm = sp_segs_of_xi(xi_eval)

   #averaging again 
   alpha_smooth_afl=0.5
   # n_avg = 5 
   for i in range(5): 
       xsm[1:-1,:] = (1.-alpha_smooth_afl)*xsm[1:-1,:] + alpha_smooth_afl*0.5*(xsm[2:,:]+xsm[:-2,:])

   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   # spacing: lower and upper back
   xsl = get_distrib(l_cut_iced[0],spac_back,spac_front,spac_int,stretch_int,force_n=force_n_lower)  # lower side (back)
   assert not ((xsl.size - 1) + (xsm.shape[0] - 1))%n_geo
   uneven = (((xsl.size - 1) + (xsm.shape[0] - 1))//n_geo) % 2
   xsu = get_distrib(l_iced[-1]-l_cut_iced[1],spac_front,spac_back,spac_int,stretch_int,force_n=force_n_upper,uneven=uneven) # upper side (back) 

   if reference_mode: 
       with h5py.File(filename_coords_base, 'w') as h5f:
           h5f.attrs.create('n_ref_upper', xsu.shape[0])
           h5f.attrs.create('n_ref_front', xsm.shape[0])
           h5f.attrs.create('n_ref_lower', xsl.shape[0])

   #join
   ICE = np.concatenate((sp(xsl[:-1]),xsm,sp(xsu[1:]+l_cut_iced[1])))

   ncor = ICE.shape[0]

   # ADDITIONAL AVERAGING: WHOLE PROFILE
   for i in range(n_avg-5): 
       ICE[1:-1,:] = (1.-alpha_smooth_afl)*ICE[1:-1,:] + alpha_smooth_afl*0.5*(ICE[2:,:]+ICE[:-2,:])

   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   # clean airfoil

   # length of polyline  
   l_clean, l_cut_clean = get_l_lr(x_afl)
   # get spline
   spa = CubicSpline(l_clean,x_afl)
   # spacing (adapt from iced)
   xsla = xsl*l_cut_clean[0]/l_cut_iced[0]
   xsua = xsu * (l_clean[-1]-l_cut_clean[1])/(l_iced[-1]-l_cut_iced[1])

   # front: linspace, then linear interpolation to l_cut
   xsmt = np.linspace(-1.,1.,n_front)
   xsma = l_cut_clean[0] + (xsmt-xsmt[0])*(l_cut_clean[1]-l_cut_clean[0])/(xsmt[-1]-xsmt[0])
   #join
   xsa = np.concatenate((xsla[:-1],xsma,xsua[1:]+l_cut_clean[1]))

   # interpolate
   AFL = spa(xsa)



   # --------------------------------------------------------------------------------
   # EXTRUDE UN-ICED STRUCTURED C-MESH 

   # blend normal extrusion (front part of aifoil) and vertical extrusion (near trailing edge) 


   XN=np.zeros((nrad,ncor,2)) #normal extrusion
   XV=np.zeros((nrad,ncor,2)) # vertical extrusion
   XN[0,:,:] = AFL
   XV[0,:,:] = AFL

   # get normal vector 
   pa = AFL[2:,:]-AFL[:-2,:]
   lp = length(pa)
   nor = np.stack((-pa[:,1] / lp, pa[:,0] / lp), axis = -1) 

   h = ncor//2 # XV: first half down, second half up

   drad = drad0

   # extrude
   for i in range(1,nrad):
       XN[i,1:-1,:] = XN[i-1,1:-1,:]+drad*nor
       XV[i,:,0] = AFL[:,0]
       XV[i,:h,1] = XV[i-1,:h,1]-drad
       XV[i,h:,1] = XV[i-1,h:,1]+drad
       drad*=stretch_fac

   a = 0.1
   o = 0.05

   # equistant spacing
   for i in range(2,nrad):
       tmp = XN[i,1:-1,:]
       ai = a*(i-1)/(nrad-2)
       tmp[:,0] = ((1.+ai)*(tmp[:,0]-o) -ai*np.abs(tmp[:,0]-o))+o
       l_old, _ = get_l_lr(tmp)
       redist = np.where(tmp[:,0]<x_max_front)[0]
       lo, up = redist[0], redist[-1]
       l_new = copy.copy(l_old)
       l_new[lo:up+1] = np.linspace(l_old[lo],l_old[up],up+1-lo)
       # Laplace smoothing for spacing
       for j in range((200*i)//nrad): 
           l_new[1:-1] = (1.-laplace_amplitude)*l_new[1:-1] + 0.5*laplace_amplitude*(l_new[2:]+l_new[:-2])
       sp = CubicSpline(l_old,tmp)
       XN[i,1:-1,:] = sp(i/nrad*l_new+(nrad-i)/nrad*l_old)

   # blend XN and XV
   X = copy.deepcopy(XN)
   t = add_dim(max_(0.,(AFL[:,0]-0.6)/0.4))
   for i in range(nrad): 
       X[i,:,:] = t*XV[i,:,:]+(1-t)*XN[i,:,:]


   # --------------------------------------------------------------------------------
   # GET d_growth FOR INCREMENTAL GROWTH 

   pa = ICE[2:,:]-ICE[:-2,:]
   lp = length(pa)
   nor = np.stack((-pa[:,1] / lp, pa[:,0] / lp), axis = -1) 

   ICE_L1 = copy.copy(ICE)
   ICE_L1[1:-1,:] += (1-t[1:-1,:])*drad0*nor
   ICE_L1[:h,1] -= t[:h,0]*drad0
   ICE_L1[h:,1] += t[h:,0]*drad0

   d_growth = np.zeros_like(X)
   d_growth[0,:,:]=(ICE-AFL)/n_steps
   d_growth[1,:,:]=(ICE_L1-X[1,:,:])/n_steps

   fac1, fac2 = np.meshgrid(d_growth[1,:,:],np.linspace(1.,0.,nrad-1))
   d_growth[1:,:,:]= (fac1*fac2).reshape(X.shape-np.array([1,0,0]))


   # --------------------------------------------------------------------------------
   # LET ICE GROW AND SMOOTH MESH 

   # functions

   # get movement factor: 1 if max_ratio is far exceeded, 0 if far below, ramp else.
   def fac_by_ramp(arr): 
       return constrain_0to1((arr-0.8*max_ratio)/(0.3*max_ratio))

   # get max dx to prevent crossing of grid lines; do not move further than this 
   # (is not very safe in realtiy...)
   def move_safely(X,dx,ldn,ldp): 
       dmax = 0.5 * min_(min_(ldp[2:-1,1:],ldp[2:-1,:-1]),min_(ldn[2:,1:-1],ldn[1:-1,1:-1]))
       ldx = length(dx)
       X[2:-1,1:-1,:] += dx*min_(1., add_dim(dmax/(ldx+1.E-16)))*fac_arr

   # plot new frame of current state
   def update_fig(): 
       fig = plt.gcf()
       fig.clf()
       ax = plt.gca()
       ax.axis('equal')
       ax.set_xticks([])
       ax.set_yticks([])
       ax.axis(ax_limits)
       for i in range(nrad):
           if not i % n_geo:
               plt.plot(X[i,:,0],X[i,:,1],'k')
           elif do_plot_inner_lines:
               plt.plot(X[i,:,0],X[i,:,1],inner_line_str+'k')
       for i in range(ncor):
           if not i % n_geo:
               plt.plot(X[:,i,0],X[:,i,1],'k')
           elif do_plot_inner_lines:
               plt.plot(X[:,i,0],X[:,i,1],inner_line_str+'k')
       fig.canvas.draw_idle()  # required to work on mpl < 1.5
       fig.canvas.flush_events()


   dradv = np.array([drad0*stretch_fac**i for i in range(nrad)])
   _, drad = np.meshgrid(np.arange(ncor),dradv)
   dradmin = drad0*(drad/drad0)**0.5

   fac = 0.5+0.5*np.cos(np.pi*constrain_0to1((ICE[:,0]-0.3)/0.7))
   fac_arr, _ = np.meshgrid(fac,np.arange(nrad))
   fac_arr = add_dim(fac_arr[2:-1,1:-1])


   # THE BIG LOOP FOR MESH MOVEMENT
   for j in range(n_steps+n_steps_after): 

       if not do_ice: 
           break

       if not j%n_steps_show and (do_show_fig or do_save_fig): 
           update_fig()

       if j < n_steps: 
           X += d_growth

       # prerequisites
       dn = X[1:,:,:]-X[:-1,:,:]
       dp = X[:,1:,:]-X[:,:-1,:]
       dp2 = X[:,2:,:]-X[:,:-2,:]
       ldp2 = length(dp2)
       dn2 = X[3:,:,:]-X[1:-2,:,:]
       ldn2 = length(dn2)
       de = np.stack((-dp2[:,:,1]*drad[:,1:-1] / ldp2, dp2[:,:,0]*drad[:,1:-1] / ldp2), axis = -1) 
       ldn = length(dn)
       ldnrat = add_dim(ldn[1:,1:-1]/ldn[:-1,1:-1])
       ldp = length(dp)
       ldprat = add_dim(ldp[2:-1,1:]/ldp[2:-1,:-1])
       an,_,_ = angle(de[:-1,:,:],dn[:,1:-1,:])

       # standard laplace
       dx = laplace_amplitude*((X[3:,1:-1,:] + stretch_fac*X[1:-2,1:-1,:])/(2.*(1.+stretch_fac)) 
               + 0.25*(X[2:-1,:-2,:] + X[2:-1,2:,:])-X[2:-1,1:-1,:])
       dx_sum = 1.*dx

       # mindist
       # normal plus (outwards)
       odn = np.stack((-dp[:,:,1] / ldp, dp[:,:,0] / ldp), axis = -1) 
       od_npp = np.sum(odn[1:-2,1:,:]*dn[1:-1,1:-1,:],axis=-1)
       od_npm = np.sum(odn[1:-2,:-1,:]*dn[1:-1,1:-1,:],axis=-1)
       dx_sum += 0.6*add_dim(max_(dradmin[2:-1,1:-1] - 0.6*od_npp,0.)**2/dradmin[2:-1,1:-1])*odn[1:-2,1:,:]
       dx_sum += 0.6*add_dim(max_(dradmin[2:-1,1:-1] - 0.6*od_npm,0.)**2/dradmin[2:-1,1:-1])*odn[1:-2,:-1,:]

       # parallel plus (right)
       odp = np.stack((dn[1:,:,1] / ldn[1:,:], -dn[1:,:,0] / ldn[1:,:]), axis = -1) 
       od_ppp = np.sum(odp[1:,:-2,:]*dp[2:-1,:-1,:],axis=-1)
       od_ppm = np.sum(odp[:-1,:-2,:]*dp[2:-1,:-1,:],axis=-1)
       dx_sum += 0.6*add_dim(max_(dpamin - 0.8*od_ppp,0.)**2/dpamin)*odp[1:,:-2,:]
       dx_sum += 0.6*add_dim(max_(dpamin - 0.8*od_ppm,0.)**2/dpamin)*odp[:-1,:-2,:]

       # parallel minus (left)
       od_pmp = np.sum(odp[1:,2:,:]*dp[2:-1,1:,:],axis=-1)
       od_pmm = np.sum(odp[:-1,2:,:]*dp[2:-1,1:,:],axis=-1)
       dx_sum += -0.6*add_dim(max_(dpamin - 0.8*od_pmp,0.)**2/dpamin)*odp[1:,2:,:]
       dx_sum += -0.6*add_dim(max_(dpamin - 0.8*od_pmm,0.)**2/dpamin)*odp[:-1,2:,:]

       # maximum width ratio
       do_pull_in = fac_by_ramp(ldnrat[:-1,:,:])
       dx = do_pull_in * (X[1:-2,1:-1,:]+max_ratio*add_dim(ldn[:-2,1:-1]/ldn2[:,1:-1])*dn2[:,1:-1,:] - X[2:-1,1:-1,:])
       dx_sum += 0.9*max_ratio_move_fac * dx

       # do_push_in = 1./ldnrat[1:,:,:] > max_ratio
       do_push_in = fac_by_ramp(1./ldnrat[1:,:,:])
       dx = do_push_in * (X[3:,1:-1,:]-dn2[:,1:-1,:]/(1+max_ratio) - X[2:-1,1:-1,:])
       dx_sum += 0.5*max_ratio_move_fac * dx

       do_push_fwd = fac_by_ramp(ldprat)
       dx = do_push_fwd * ((max_ratio*X[2:-1,:-2,:]+X[2:-1,2:,:])/(1.+max_ratio) - X[2:-1,1:-1,:])
       dx_sum += 0.5*max_ratio_move_fac * dx

       do_push_back = fac_by_ramp(1./ldprat)
       dx = do_push_back * ((max_ratio*X[2:-1,2:,:]+X[2:-1,:-2,:])/(1.+max_ratio) - X[2:-1,1:-1,:])
       dx_sum += 0.5*max_ratio_move_fac * dx

       # angle
       xefac = (0.005+0.12*add_dim(constrain_0to1((an[1:-1,:]-0.9*alpha_max)/(0.2*alpha_max))))
       dx = xefac*(de[1:-2,:,:] - dn[1:-1,1:-1])
       dx_sum += dx

       move_safely(X,0.5*dx_sum,ldn,ldp)

   # --------------------------------------------------------------------------------
   # RESAMPLE TO NEW GRID

   post_dradv = np.array([post_drad0*post_stretch**i for i in range(post_nrad)])
   post_dradv = min_(post_dradv,post_dradmax)

   radv = np.concatenate(([0.,],np.cumsum(dradv[:nrad-1])))
   post_radv = np.concatenate(([0.,],np.cumsum(post_dradv[:post_nrad-1])))
   # Post mesh must be thinner than original mesh
   assert post_radv[-1] <= radv[-1] 





   # normal
   dn = X[1:,:,:]-X[:-1,:,:]
   ldn = length(dn)
   xi_rad = np.concatenate((np.zeros_like(X[0:1,:,0]),np.cumsum(ldn,axis=0)),axis=0)
   xi_rad_eval = np.ndarray((post_nrad,ncor))
   for i in range(ncor): 
       xirad_inter = interp1d(radv,xi_rad[:,i])
       xi_rad_eval[:,i] = xirad_inter(post_radv)
   for i in range(1,n_geo): 
       xi_rad_eval[i::n_geo,:] = (1.-np.float(i)/n_geo) * xi_rad_eval[0:-1:n_geo,:] + (np.float(i)/n_geo)*xi_rad_eval[n_geo::n_geo,:]
   XTMP = X
   X = np.ndarray((post_nrad,ncor,2))
   for i in range(ncor): 
       x_inter = interp1d(xi_rad[:,i],XTMP[:,i,:],axis=0)
       X[:,i,:] = x_inter(xi_rad_eval[:,i])

   # adjust indices
   n_elems_rad = post_n_elems_rad
   nrad = post_nrad
   radv = post_radv

   # parallel
   dp = X[:,1:,:]-X[:,:-1,:]
   ldp = length(dp)
   xi_cor = np.concatenate((np.zeros_like(X[:,0:1,0]),np.cumsum(ldp,axis=1)),axis=1)
   xi_cor_eval = copy.copy(xi_cor)
   for i in range(1,n_geo): 
       xi_cor_eval[:,i::n_geo] = (1.-np.float(i)/n_geo) * xi_cor_eval[:,0:-1:n_geo] + (np.float(i)/n_geo)*xi_cor_eval[:,n_geo::n_geo]
   for i in range(nrad): 
       x_inter = interp1d(xi_cor[i,:],X[i,:,:],axis=0)
       X[i,:,:] = x_inter(xi_cor_eval[i,:])


   # --------------------------------------------------------------------------------
   # READ REF COORDS AND MOVE TO REF COORDS
   if reference_mode: 
       with h5py.File(filename_coords_base, 'r+') as h5f:
           h5f.create_dataset('coords', data=X)
   else:
       with h5py.File(filename_coords_base, 'r') as h5f:
           X_ref = np.array(h5f['coords'])

       start_blend = 0.2
       end_blend = 1.0
       for i in range(nrad): 
           # linblend = constrain_0to1(((i+1.)/nrad-start_blend)/(end_blend-start_blend))
           linblend = constrain_0to1((radv[i]/radv[-1]-start_blend)/(end_blend-start_blend))
           # fac = 0.5-0.5*np.cos(np.pi*linblend**2)
           # fac = -2*linblend**3+3*linblend**2
           fac = linblend**2.
           X[i,:,:] = (1.-fac)*X[i,:,:] + fac*X_ref[i,:,:]


   if do_show_fig or do_save_fig: 
       update_fig()

   if do_show_fig:
       plt.ioff()
       plt.show()

   if do_save_fig: 
       outdir = "fig_output/meshes/global_laplace_"
       i = len(glob.glob(outdir+"*"))
       plt.savefig(outdir+str(i+1)+".png")


   # --------------------------------------------------------------------------------
   # BUILD STRUCTURED WAKE MESH 

   al_upper = al_wake + phi_wake
   al_lower = al_wake - phi_wake

   b_upper = np.array(circ_ctr[0:2])+r_circ*np.array([np.cos(al_upper), np.sin(al_upper)])
   b_lower = np.array(circ_ctr[0:2])+r_circ*np.array([np.cos(al_lower), np.sin(al_lower)])
   b_ctr = (b_upper+b_lower)/2.


   sum_wake = spac_back
   n_wake = 1
   while (sum_wake < r_circ) or (n_wake % n_geo): 
       sum_wake += spac_back*stretch_fac_wake**n_wake
       n_wake +=1 

   wakev = np.cumsum(np.array([(stretch_fac_wake**i) for i in range(n_wake)]))*spac_back/sum_wake
   radv /=radv[-1]

   def blow(v1,v2): 
       return np.matmul(np.transpose(np.array([v1])),np.array([v2]))
   XL = np.ndarray((nrad,n_wake,2))
   XL[0,:,:] = blow((1.-wakev),X[0,0,:]) + blow(wakev,b_ctr)
   XL[-1,:,:] = blow((1.-wakev),X[-1,0,:]) + blow(wakev,b_lower)
   for idim in range(2): 
       XL[:,:,idim] = blow((1.-radv),XL[0,:,idim]) + blow(radv,XL[-1,:,idim])

   XU = np.ndarray((nrad,n_wake,2))
   XU[0,:,:] = blow((1.-wakev),X[0,-1,:]) + blow(wakev,b_ctr)
   XU[-1,:,:] = blow((1.-wakev),X[-1,-1,:]) + blow(wakev,b_upper)
   for idim in range(2): 
       XU[:,:,idim] = blow((1.-radv),XU[0,:,idim]) + blow(radv,XU[-1,:,idim])

   XE = np.concatenate((np.flip(XL,axis=1),X,XU),axis = 1)

   # smooth transition to wake (prevents triangles in gmsh)
   halfwidth = 10*n_geo
   depth = int(round(post_nrad*0.8-0.5))
   for i in range(1,depth+1): 
       alpha = 0.5*((depth+1-i)/depth)
       for hw in range(1,halfwidth + 1): 
           for c in [n_wake, XE.shape[1]-n_wake]:
               XE[-i,c-hw:c+hw+1,:] = (1.-alpha)*XE[-i,c-hw:c+hw+1,:] + alpha*0.5*(XE[-i,c-hw-1:c+hw,:]+XE[-i,c-hw+1:c+hw+2,:])


   assert XE.shape[1] % (2*n_geo) == 1
   for i in range(1,n_geo): 
       XE[-1,i::n_geo,:] = (i*XE[-1,n_geo::n_geo,:]+(n_geo-i)*XE[-1,:-n_geo:n_geo])/n_geo

   # --------------------------------------------------------------------------------

   if not do_write_mesh: 
       sys.exit()


   if reference_mode: # build unstructured gmsh

       # GET CIRCULAR OUTER UNSTRUCTURED GMSH 
       geom = pygmsh.built_in.Geometry()
       geom.add_raw_code("Mesh.Algorithm = 9;")

       XM = np.concatenate((XE[-1,::n_geo,:],0.*XE[-1,::n_geo,0:1]),axis=-1)

       def to3d(v): 
           return np.concatenate((v,np.array([0.])))

       n_circ += (XM.shape[0]+n_circ)%2
       for phi in np.linspace(al_upper,al_lower+2.*np.pi,n_circ+2)[1:-1]:
           XM = np.concatenate((XM,circ_ctr + r_circ * np.array([[np.cos(phi),np.sin(phi),0.]])),axis=0)

       poly = geom.add_polygon(XM)


       geom.set_transfinite_lines(poly.lines, 1)

       geom.add_raw_code("Recombine Surface {{{}}};".format(poly.surface.id))


       zunit = [0., 0., 1.]
       top,vol,lat = geom.extrude(poly.surface,translation_axis=zunit,num_layers=1,recombine=True)

       geom.add_physical(vol)
       # geom.add_physical(poly.surface,label="BOTTOM")
       # geom.add_physical(top,label="TOP")
       # geom.add_physical(lat[:-n_circ-1],label="INNER")
       geom.add_physical(lat[-n_circ-1:],label="BC_OUTER")

       geom.add_raw_code("Mesh.MshFileVersion = 2.2;")
       mesh = pygmsh.generate_mesh(geom,geo_filename=filename_geo)#,msh_filename=filename_msh)
       if "wedge" in mesh.cells: 
           print("\nERROR! {} triangles created!!!\n".format(len(mesh.cells["wedge"])))
           sys.exit()

       meshio.write(filename_msh,mesh,file_format='gmsh2-ascii')
       print("msh file: "+filename_msh)
       print()
       print()

   else: # build FLEXI mesh

       # write coords
       with h5py.File(filename_coords, 'w') as h5f:
            h5f.create_dataset('coords', data=XE)
            h5f.attrs.create("n_wake", data=n_wake//n_geo)

       # copy and complete hopr ini file
       with open(hoprbasefile,'r') as pf: 
           prms = pf.read()
       prms += "\nprojectname = "   +name_str
       prms += "\nfilename = "      +filename_msh
       prms += "\nCoordsFilename = "+filename_coords
       hoprfile = name_str+"_parameter_hopr.ini"
       with open(hoprfile,'w') as pf: 
           pf.write(prms)

       # execute hopr
       print("Run hopr...")
       args = [hopr_path, hoprfile]
       output=subprocess.run(args,stdout=subprocess.PIPE)
       stdout_str=output.stdout.decode("utf-8")
       print("Run hopr DONE")
       return stdout_str

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

if __name__ == "__main__": 
   hopr_path   = sys.argv[1]
   # reference mode: save reference coordinates and build unstructured mesh 
   # normal mode: bend to ref coords and execute hopr with reference unstructured mesh
   reference_mode = (hopr_path == "NONE")
   hoprbasefile = sys.argv[2]
   name_str    = sys.argv[3]
   random_vars = sys.argv[4].split("_") # random distribution; normal(0.,1.,1) is correct, others improve stability
   n_avg = int(sys.argv[5])

   hopr_stdout_str = make_mesh(hopr_path,reference_mode,hoprbasefile,name_str,random_vars,n_avg)
   if hopr_stdout_str: 
       print(hopr_stdout_str)
