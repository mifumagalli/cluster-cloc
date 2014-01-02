from parse import sfrl1_parse, mv2ln
from fall_cpp_wrapper import fall_cpp_wrapper
import os
from numpy import *
from scipy import *
from scipy.stats import kstest, uniform
import struct
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from scipy.ndimage.filters import gaussian_filter
from operator import itemgetter
from scipy.optimize import curve_fit, leastsq, fmin_tnc
from mpfit import *
import numpy.oldnumeric as Numeric
from datetime import datetime
import time



rc('text', usetex=True)
rc('font', family='serif')
data =( (102.3, -13.8),\
        (32.2,-16.3),\
        (4.8,-14.7),\
        (15.0,-14.9),\
        (13.3,-13.3),\
        (57.6,-14.5),\
        (58.1, -15.6),\
        (1.3, -11.4),\
        (0.46, -11.8),\
        (1.0, -11.4),\
        (5.2, -12.8),\
        (2.1, -12.3),\
        (2.6, -12.9),\
        (0.16,-10.9),\
        (0.14, -9.7),\
        (0.24, -11.0),\
        (0.52, -10.1),\
        (0.63, -11.9),\
        (2.7, -13.0),\
        (0.78, -11.8),\
        (0.18, -9.1),\
        (0.077, -9.7),\
        (0.2,-11.3),\
        (0.78, -11.2),\
        (0.31, -10.9),\
        (0.054,-14.0))
log_sfr = array([log10(data[i][0]) for i in xrange(len(data))] )
lnMV = mv2ln(array([data[i][1] for i in xrange(len(data))] ))

out, table = sfrl1_parse('grid.dat')

#data_xind = array([min(enumerate(abs(log_sfr[i] - out.sfr_x)), key=itemgetter(1))[0] \
#                   for i in xrange(len(log_sfr)])
#data_yind = array([min(enumerate(abs(lnMV[i] - out.x)), key=itemgetter(1))[0] \
#                   for i in xrange(len(MV)])


data_xind = array( [argmin(abs(log_sfr[i] - out.sfr_x)) for i in xrange(len(log_sfr))])

#data_yind = array( [argmin(abs(lnMV[i] - out.x)) for i in xrange(len(log_sfr))])
data_yind = lnMV

def data2q(dat_sfr_ind, dat_mv_ind, tab_out, tab_tab):
    '''use data and model to compute quantile plot'''
    q = zeros(len(dat_mv_ind))
    for i in xrange(len(dat_mv_ind)):
        pdf = tab_tab[dat_sfr_ind[i]]   
        csum = cumsum(pdf)
        csum /= csum[-1]
        #q[i] = csum[dat_mv_ind[i]]
        q[i] = interp(dat_mv_ind[i],tab_out.x , csum)
    return q


def call_wrap2(p, mmax=1e6, fc=0.05, x=None, y=None, err=None, fjac=None, more_out=False):
        #print mmax
        #print fc
        print p
        #print shape(p)
        nfake = 10

        fall_cpp_wrapper(t0 = p[0], t1=p[0]+p[1], mmax = mmax, f_c=fc,\
           mmin=10.**p[2], cmfslope = p[3], age_slope = p[4], \
           grid_out = 'bast.dat')

        out, table = sfrl1_parse('bast.dat')
      
        qq = data2q(data_xind, data_yind, out, table)
        #print qq
        D, p = kstest(qq, 'uniform')
        print [1-p, D]
        print '------'
        status = 0

        #return [status, (1-p)*ones(nfake)]
        if more_out == True:
            var_rat = std(qq)**2*12
            mean_dat = mean(qq)
            return [D, p, var_rat, mean_dat]
        return [status, D*ones(nfake)]

def fit_mod(fc=0.01, mmax = 1e7):
    t0 = 0.1
    delta_t = 0.5
    mmin = log10(200.)
    age_slope = -0.9
    cmf_slope = -2.0
    starting = array([t0, delta_t, mmin, cmf_slope,  age_slope])
    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.], \
                'step':0.} 
                                                for i in range(5)]
    fac = 0.5
    print starting
    parinfo[0]['limited'][0]=1
    parinfo[0]['limited'][1]=1
    parinfo[0]['limits'][0]=0.
    parinfo[0]['limits'][1]=10.
    parinfo[0]['value']=0.1
    parinfo[0]['step']=0.01*fac

    parinfo[1]['limited'][0]=1
    parinfo[1]['limited'][1]=1
    parinfo[1]['limits'][0]=0.
    parinfo[1]['limits'][1]=10.
    parinfo[1]['value']=0.5
    parinfo[1]['step']=0.01*fac

    parinfo[2]['limited'][0]=1
    parinfo[2]['limited'][1]=1
    parinfo[2]['limits'][0]=2.
    parinfo[2]['limits'][1]=4.
    parinfo[2]['value']=2.2
    parinfo[2]['step']=0.01*fac
    
    parinfo[3]['limited'][0]=1
    parinfo[3]['limited'][1]=1
    parinfo[3]['limits'][0]=-2.2
    parinfo[3]['limits'][1]=-1.8
    parinfo[3]['value']=-2.
    parinfo[3]['step']=0.01*fac
    
    parinfo[4]['limited'][0]=1
    parinfo[4]['limited'][1]=1
    parinfo[4]['limits'][0]=-0.2
    parinfo[4]['limits'][1]=0.2
    parinfo[4]['value']=0
    parinfo[4]['step']=0.01*fac
    
    sss = call_wrap2(starting)
    return mpfit(call_wrap2, parinfo = parinfo, \
                    functkw={'mmax':mmax, 'fc':fc})
    #return curve_fit(call_wrap_maker(fc, mmax) , zeros(10.), zeros(10.), starting )
    #return leastsq(call_wrap_maker(fc, mmax), starting , \
    #               diag = array([1,1,100,1,1])*1000, epsfcn = 0.1, ftol=0.0005)
    #return fmin_tnc(call_wrap_maker(fc, mmax), starting, approx_grad=True, \
     #  bounds=( (0., 5.), (0.,5.), (10., 10000.), (-2.5, -1.5), (-1.2, 0.2)))


logfc = linspace(-3., 0., 20)
logmmax = linspace(5, 9, 20)
nmodels = len(logfc)*len(logmmax)

ff = str(long(nmodels))+'f8'
output = recarray((1,),dtype=[('logmmax', ff), ('logfc',ff), ('t0',ff),\
                     ('delta_t',ff), ('mmin',ff), \
                     ('cmf_slope',ff), ('age_slope', ff), ('p',ff),('D', ff), ('var_rat',ff),\
                     ('mean_dat', ff)])[0]

name = 'bast.npy'
k = 0
for lfc in logfc:
    for lmmax in logmmax:
        output.logmmax[k] = lmmax
        output.logfc[k] = lfc
        output.t0[k] = 0
        output.delta_t[k] = 0
        output.mmin[k] = 0
        output.cmf_slope[k] = 0
        output.age_slope[k] = 0
        output.p[k] = 0
        output.D[k] = 0
        output.var_rat[k] = 0
        output.mean_dat[k] = 0
        k = k +1
save(name, output)
k = 0
for lfc in logfc:
    for lmmax in logmmax:
        p0 = fit_mod(fc = 10.**lfc, mmax = 10.**lmmax)
        D, p, var_rat, mean_dat = call_wrap2(p0.params, mmax = 10.**lmmax, \
                                              fc=10.**lfc, more_out=True)
        output.logmmax[k] = lmmax
        output.logfc[k] = lfc
        output.t0[k] = p0.params[0]
        output.delta_t[k] = p0.params[1]
        output.mmin[k] = p0.params[2]
        output.cmf_slope[k] = p0.params[3]
        output.age_slope[k] = p0.params[4]
        output.p[k] = p
        output.D[k] = D
        output.var_rat[k] = var_rat
        output.mean_dat[k] = mean_dat
        k = k + 1
        save(name, output)

save(name, output)

#y0 = 1-call_wrap_maker(fc=0.01, mmax=1e7)(p0[0])[0]
#fall_cpp_wrapper()
#qq = data2q(data_xind, data_yind, out, table)
#var_rat = std(qq)**2*12
#mean_dat = mean(qq)
#D, p = kstest(qq, 'uniform')



