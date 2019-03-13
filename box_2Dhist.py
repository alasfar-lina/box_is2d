import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import uproot
# from pyik.numpyext import centers

def hist2D_norm (arr, bins):
    return   np.sum(arr*np.diff(bins[-1].T))

def make_box (i,j,narrxy,bins,vmin,vmax,scale):
    normx = np.sum(narrxy,axis=1)
    arrxy= narrxy/normx
    x= bins[0]
    y=bins[1]
    dx= np.diff(x)
    dy=np.diff(y)
    narrx,narry = (0,0)
    norm= hist2D_norm(arrxy, bins)
    if scale == 'linear':
        if vmin == None and vmax == None:
            narrx= np.sqrt( ( arrxy[i,j]  )  / norm  ) * dx[i]
            narry= np.sqrt( ( arrxy[i,j]  )  / norm  ) * dy[j]
        else:
            if arrxy[i,j] > vmax :
                arrxy[i,j] = vmax
            if arrxy[i,j] < vmin :
                arrxy[i,j]= vmin

            narrx= ( ( arrxy[i,j] - vmin )  / vmax ) * dx[i]
            narry= ( ( arrxy[i,j] - vmin )  / vmax ) * dy[j]
    if scale == 'log':
        # if vmin == None and vmax == None:
        #     narrx= (1-0.5*( np.log10( arrxy[i,j]/norm  ) ))* dx[i]
        #     narry= (1-0.5*( np.log10( arrxy[i,j]/norm  )) )* dy[j]
        # else:
        if arrxy[i,j] > vmax :
            arrxy[i,j] = vmax
        if arrxy[i,j] < vmin:
            arrxy[i,j]= vmin
        # narrx = (dx[i] *( np.log10(arrxy[i,j]/vmax) ))
        # narry = dy[i] *( np.log10(arrxy[i,j]/vmax) )

        # narrx= 0.5*dx[i]-( np.log10(arrxy[i,j] / vmin)) *0.5*dx[i]
        # narry= 0.5*dy[j]-(np.log10(arrxy[i,j] /vmin) ) ) *0.5*dy[j]
        narrx= (1- np.log( arrxy[i,j]  )  / np.log(vmin) ) * dx[i]
        narry= (1- np.log( arrxy[i,j]  )  / np.log(vmin) ) * dy[j]
    return Rectangle((x[i] +(dx[i]-(narrx))/2.0 , y[j]+(dy[j]-(narry))/2.0 ),(narrx),(narry),fill=False)

def root_box_plot(ax, x_bins, y_bins, arr2D, vmin=None, vmax=None, scale ='linear',**kwargs):
    """
    Plots a ROOT 2D histogram in box get_style
    the area of the box is equal to the value of the bin/norm of the histogram.

    Parameters
    ----------
    ax: plot axis
    matplotlib.pyplot.axis
    x_bins : x bin edges
    y_bins, y bin edges
    arr2D: bin contents
    vmin : minimum value, zero size box
    vmax : maximum value, full bin width box
    **kwagrs (facecolor , edgecolor alpha... )

    Returns
    -------
    PatchCollection of rectagles:
    a box 2D histogram

    Dependencies
    --------
    numpy, matplotlib, uproot.
    if you do not have root histograms, use Hist2D from
    rootpy.plotting to create one
    """

    boxes =[]
    m = x_bins.shape[-1]-1
    n=y_bins.shape[-1]-1
    for j in range(n): # over x-axis
        for i in range(m): #over y-axis
            box=make_box(i,j,arr2D,bins,vmin,vmax,scale)
            boxes.append(box)
    pc = PatchCollection(boxes ,**kwargs)
    ax.add_collection(pc)
    plt.xlim(bins[0][0]-0.1,bins[0][-1])
    plt.ylim(bins[1][0]-0.1,bins[1][-1])


# Example 1
#
# from rootpy.plotting import Hist2D
# from root_numpy import array2hist
# hist_e = Hist2D(6, 2,5 , 16, 0, 10, type='D')
# arr_e = np.loadtxt("/home/lina/charged-hadron-spectra/data/eff.dat")
#
# histo_epsilon= array2hist(arr_e, hist_e)
#
# ax = plt.axes()
# root_box_plot(ax,histo_epsilon, vmin = 0.1, vmax=2 ,scale="linear")
#
# plt.xlabel(r"$\eta$")
# plt.ylabel(r"$p_T$ [GeV/c]")
# plt.title("efficiency box histogram ")
# plt.show()
#Example 2
# from rootpy.io import root_open
# from rootpy.plotting.style import get_style, set_style
# import logging
# mpl_logger = logging.getLogger('matplotlib')
# mpl_logger.setLevel(logging.WARNING)
# set_style('lhcb', mpl=True)
#
# # Example 2
#
f_rmu= uproot.open("/Users/lina/charged-hadron-spectra/bin_to_bin/root_files/response_matrix_magdown_MC.root")
histo_rm1u= f_rmu['eta_b1']  #2,2.5
arr, bins  = (histo_rm1u).numpy()
x_bins= bins[0]
y_bins= bins[1]

m = x_bins.shape[-1]-1
n=y_bins.shape[-1]-1
normx = np.sum(arr,axis=1)
arr= arr/normx
for j in range(n): # over x-axis
    for i in range(m): #over y-axis
        if arr[i,j] > 1:
            arr[i,j] = 1
        if arr[i,j] <0.01:
            arr[i,j] =0.01
ax = plt.axes()
root_box_plot(ax,x_bins, y_bins, arr, vmin=1e-2, vmax= 2,scale='log'  ,edgecolor='crimson',facecolor='slateblue',alpha=0.5)

ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
plt.xlabel(r" Truth  matched $p_T$ [GeV/c] ")
plt.ylabel(r"Rec $p_T$ [GeV/c]")
plt.title(r"Response matrix , $\eta \in [2.0,\,2.5)$  ")
plt.show()
