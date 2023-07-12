#####################################################################################
# 
# Helper functions used by other scripts
# 
#####################################################################################

from scipy.interpolate import CubicSpline
import numpy as np
import csv 

def interpolate(arr,spac=None,nseg=None,spminmax=None,xs_in=None,getpa=False,getnor=False,getxs=False): 
    l = np.cumsum( np.linalg.norm(arr[1:,:]-arr[:-1,:], axis = 1))
    l = np.concatenate((np.array([0]),l)) 
    l, ui = np.unique(l,return_index=True)
    arr = arr[ui,:]
    sp = CubicSpline(l,arr)
    if spac: 
        nseg = int(l[-1]/spac)+1
        xs = np.linspace(0.,l[-1],nseg+1)
    elif nseg: 
        xs = np.linspace(0.,l[-1],nseg+1)
    elif xs_in is not None: 
        xs=xs_in
    arr_out = sp(xs)
    returnargs = [arr_out]
    if getpa or getnor: 
        pa = sp(xs, nu=1)
        lp = np.linalg.norm(pa, axis = 1)
    if getpa: 
        returnargs.append( np.stack((pa[:,0] / lp, pa[:,1] / lp), axis = -1) )
    if getnor: 
        returnargs.append( np.stack((-pa[:,1] / lp, pa[:,0] / lp), axis = -1) )
    if getxs: 
        returnargs.append( xs )
    if len(returnargs) == 1: 
        return returnargs[0]
    else:
        return returnargs

def linlen(arr): 
    l = np.cumsum( np.linalg.norm(arr[1:,:]-arr[:-1,:], axis = 1))
    return np.concatenate((np.array([0]),l)) 


def read_coords(path):
    with open(path) as csvfile:
        reader = csv.reader(csvfile)
        return np.array([[float(e) for e in r] for r in reader])

# --------------------------------------------------------------------------------

def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def do_intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def x_intersect(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]
    div = det(xdiff, ydiff)
    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return np.array([[x, y]])

def mindist(pt, array): 
    return np.argmin((pt[0]-array[:,0])**2 + (pt[1]-array[:,1])**2)

# --------------------------------------------------------------------------------

def path_union(x_afl,x_ice): 
    n_ice = x_ice.shape[0]
    x_out = np.ndarray((0,2))

    # search nearest n(i) for each point i 
    j_mindist = np.arange(n_ice)
    for i in range(n_ice): 
        j_mindist[i] = mindist(x_ice[i,:],x_afl)

    j_min = np.minimum(j_mindist[1:],j_mindist[:-1])-1
    j_max = np.maximum(j_mindist[1:],j_mindist[:-1])+1

    # check intersection for i:i+1 and min(n(i),n(i+1))-1 : max(n(i),n(i+1)+1
    # change sign for each intersection
    inside = True
    j_start = 0
    segs = []
    for i in range(n_ice-1): 
        for j in range(j_min[i],j_max[i]):
            if do_intersect(x_ice[i,:],x_ice[i+1,:],x_afl[j,:],x_afl[j+1,:]): 
                if inside: 
                    x_out=np.concatenate((x_out,x_afl[j_start:j+1,:]))
                    i_start = i+1
                    inside = False
                else: 
                    x_out=np.concatenate((x_out,x_ice[i_start:i+1,:]))
                    j_start   = j+1
                    inside = True
                x_out=np.concatenate((x_out,x_intersect(x_ice[i:i+2,:],x_afl[j:j+2,:])))
    x_out=np.concatenate((x_out,x_afl[j_start:,:]))
    return x_out
