import os
from numpy import *
from scipy import *
import struct
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from scipy.ndimage.filters import gaussian_filter


def fall_cpp_wrapper( t0=0.0,t1=0.1, mmin=100, mmax=1e9,\
        cmfslope=-2.0, gamma_min=5.05e18, gamma_max=1.198e20,\
        eta=-0.688, age_slope=-0.9, start=15., stop=35.,\
        step=0.125, sfr_err=0.5, sfr_start=-4, sfr_stop=4,\
        length=1e9, f_c=0.01, grid_out='grid.dat', obs_err=0.4604, sfr_step = 0.0125):
    s=[t0, t1, mmin, mmax, cmfslope, gamma_min, gamma_max, eta, age_slope,\
       start, stop, step, length, f_c, obs_err, sfr_start, sfr_stop,\
       sfr_step,sfr_err]
    st = ''
    #print s
    for si in s:
        #print si
        st = st + ' ' + str(si)
    dir_ = '/Users/rdasilva/Dropbox/sfr_l1_plots/calc/'
    command = dir_+'fall '+st+' '+grid_out
    print command
    os.system(command)


