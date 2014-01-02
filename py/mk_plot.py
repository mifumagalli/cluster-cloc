import os
from numpy import *
from scipy import *
import struct
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import rc
import pdb
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from scipy.ndimage.filters import gaussian_filter

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}

rc('font', **font)
rc('text', usetex=True)
rc('font', family='serif')

    #file00='dat/grid_def.dat'
    #file='dat/grid.dat'
def mk_plot(img, xr, yr, title, figname, irange=[None,None]):
    fig = plt.figure(figsize=(9,6))
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel(r'$\log_{10} {\cal F}_c $')
    ax1.set_ylabel(r'$\log_{10} m_{\rm max} \,  [M_\odot]$')
    ax1.set_title(title)
    img = ax1.imshow(img, extent=[xr[0], xr[1], yr[0], yr[1]], aspect='auto',\
              vmin=irange[0], vmax=irange[1])
    img.set_cmap('spectral')
    cb = plt.colorbar(img)
    for l in cb.ax.yaxis.get_ticklabels():
        l.set_family("Normal")
    #plt.colorbar()
    #ax1.annotate(tt,[0.05,0.85],xycoords='axes fraction', ha='left', \
    #    size = 'medium' )
    #ax1.set_xlim([-30, 0])
    #plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
    fig.savefig(figname)
    plt.show()

def plotwrap(prefix, file):
    o = load(file)
    xs = unique(o['logfc'])
    ys = unique(o['logmmax'])
    
    xr = [min(xs), max(xs)]
    yr = [min(ys), max(ys)]
    
    p_img = reshape(o['p'], [len(xs), len(ys)])[::-1,:]
    D_img = reshape(o['D'], [len(xs), len(ys)])[::-1,:]
    meandat_img = reshape(o['mean_dat'], [len(xs), len(ys)])[::-1,:]
    varrat_img = reshape(o['var_rat'], [len(xs), len(ys)])[::-1,:]
    mk_plot(p_img, xr, yr, r'K-S $p$-value', prefix+'ks.eps', \
           irange=[0,0.7116709])
    mk_plot(D_img, xr, yr, r'K-S $D$ Statistic', prefix+'D.eps',\
           irange=[0,1.])
    mk_plot(meandat_img, xr, yr, r'Mean of Quantiles', prefix+'mean_dat.eps',\
          irange=[0.,1.])
    mk_plot(varrat_img, xr, yr, r'Scaled Variance of Quantiles', \
                                   prefix+'var_rat.eps', irange=[0,2.5])

prefix = ''
plotwrap(prefix, 'grid.npy')

prefix = 'bast_'
plotwrap(prefix, 'bast.npy')
