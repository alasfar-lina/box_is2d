import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from root_numpy import hist2array
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
    if scale == "linear":
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
    if scale == "log":
        if vmin == None and vmax == None:
            narrx= (1-0.5*( np.log10( arrxy[i,j]/norm  ) ))* dx[i]
            narry= (1-0.5*( np.log10( arrxy[i,j]/norm  )) )* dy[j]
        else:
            if arrxy[i,j] > vmax :
                arrxy[i,j] = vmax
            if arrxy[i,j] < vmin:
                arrxy[i,j]= vmin
        narrx= (1 - np.log(arrxy[i,j]) / np.log(vmin)) *dx[i]
        narry= (1- np.log(arrxy[i,j]) / np.log(vmin) ) *dy[j]
            # narrx= ( np.log( arrxy[i,j]  )  / np.log(vmax) ) * dx[i]
            # narry= ( np.log( arrxy[i,j]  )  / np.log(vmax) ) * dy[j]

    return Rectangle((x[i],y[j]),(narrx),(narry),fill=False)

def root_box_plot(ax,hist2D,vmin=None,vmax=None, scale ='linear',**kwargs):
    """
    Plots a ROOT 2D histogram in box get_style
    the area of the box is equal to the value of the bin/norm of the histogram.

    Parameters
    ----------
    ax: plot axis
    matplotlib.pyplot.axis
    hist2D: ROOT TH2 object (2 histogram)
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
    arr2D, bins = hist2array(hist2D,return_edges= True)  # hist2D.numpy()
    x_bins= bins[0]
    y_bins= bins[1]
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
# f_rmu= root_open("~/charged-hadron-spectra/root_files/response_matrix_magdown_MC.root")
# histo_rm1u= f_rmu.eta_b1  #2,2.5
# ax = plt.axes()
# root_box_plot(ax,histo_rm1u, vmin=0.01, vmax= 10 , scale= 'log' ,facecolor='C4',edgecolor='C3',alpha=0.5)
# arr, bins  = hist2array(histo_rm1u, return_edges= True)
#
#
# x_bins= bins[0]
# y_bins= bins[1]
#
# m = x_bins.shape[-1]-1
# n=y_bins.shape[-1]-1
# normx = np.sum(arr,axis=1)
# arr= arr/normx
# for j in range(n): # over x-axis
#     for i in range(m): #over y-axis
#         if arr[i,j] > 1:
#             arr[i,j] = 1
#         if arr[i,j] <0.01:
#             arr[i,j] =0.01
#
#
# plt.xlabel(r" Truth  matched $p_T$ [GeV/c] ")
# plt.ylabel(r"Rec $p_T$ [GeV/c]")
# plt.title(r"Response matrix , $\eta \in [2.0,\,2.5)$  ")
# plt.show()
