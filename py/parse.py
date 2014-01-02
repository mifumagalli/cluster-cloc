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

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}

rc('font', **font)
def ln2mv(x):
    oo = -2.5*log10(exp(x)/1.1967453e40)-48.6
    return oo
def mv2ln(x):
    oo = log( 10.**((x +48.6)/(-2.5))*1.1967453e40)
    return oo
def ln2mvpdf(x):
    f = abs(-2.5*log10(exp(5))/5)
    return x/f

def read_all_doubles(name):
    ''' read in a file of all binary doubles '''
    # start by getting the size of the file in bytes
    size = os.stat(name).st_size
    if (size % 8) !=0:
        print "Size Did not match a multiple of 8!!!!"
    # open the file, read it, and put all the values in the right place
    vals = struct.unpack('d'*(size/8),open(name).read(size))    
    return array(vals)

def sfrl1_parse(file):
    doubs = read_all_doubles(file)
    cursor = 0

    sfr_len = doubs[0]
    cursor = 1

    sfr_x = doubs[cursor: sfr_len + cursor]
    cursor += sfr_len 

    x_len = doubs[cursor]
    cursor += 1

    x = doubs[cursor: x_len + cursor]
    cursor += x_len

    pdf_l_obs = doubs[cursor: x_len + cursor]
    cursor += x_len

    pdf_l = doubs[cursor: x_len + cursor]
    cursor += x_len

    table = doubs[cursor:]
    table = reshape(table, (sfr_len, x_len))
    f = 'f8'
    sxl = str(long(x_len))
    sfrl = str(long(sfr_len))
    ff=sxl+f
    fff = sfrl + f
    out = recarray((1,),dtype=[('x',ff),('pdf_l',ff), ('pdf_l_obs',ff), \
                     ('sfr_x',fff), ('q95', fff), ('q50',fff),('q5', fff)])[0]
    out.x[:] = x[:]
    out.pdf_l[:] = pdf_l[:]
    out.pdf_l_obs[:] = pdf_l_obs[:]
    out.sfr_x[:] = sfr_x[:]
    q95 = zeros(sfr_len)
    q50 = zeros(sfr_len)
    q5 = zeros(sfr_len)
    for i in arange(sfr_len):
        this_pdf = table[i]
        csum = cumsum(this_pdf)
        csum /= csum[-1]
        q95[i] = x[argmin(abs(0.95-csum))]
        q50[i] = x[argmin(abs(0.5-csum))]
        q5[i] = x[argmin(abs(0.05-csum))]        
    out.q95[:] = q95[:]
    out.q50[:] = q50[:]
    out.q5[:] = q5[:]
    #pdb.set_trace()
    #print [out.x[0:3], x[0:3]]
    return out, table

if __name__  == '__main__':
    rc('text', usetex=True)
    rc('font', family='serif')

    #file00='dat/grid_def.dat'
    #file='dat/grid.dat'
    def mk_plot(file00, file, tt, this_err):
        global k
        global def_err
        name = file[0:-4]
        out, table = sfrl1_parse(file)
        out00, table00 = sfrl1_parse(file00)
        
        pdf00=out00.pdf_l
        
        
        alpha = 0.3
        ax1 = fig.add_subplot(gs1[k,0])
        #ax1.fill_between(ln2mv(out.x), \
        #               ln2mvpdf(out.pdf_l), 0, color=appa('#FFFFFF',alpha))

        convolved_pdf_l = gaussian_filter(out.pdf_l, this_err/(out00.x[1]-out00.x[0]))
        convolved_pdf_l00 = gaussian_filter(out00.pdf_l, def_err/(out00.x[1]-out00.x[0]))
        ax1.fill_between(ln2mv(out00.x), \
                       ln2mvpdf(convolved_pdf_l00), 0, color=appa('#FF0000',alpha))
        ax1.plot(ln2mv(out.x), ln2mvpdf(out.pdf_l), color=appa('#000000', 1.0))
        ax1.plot(ln2mv(out.x), ln2mvpdf(convolved_pdf_l), color='red')

        ax1.set_xlabel(r'$L \, [M_V]$')
        ax1.set_ylabel(r'$p(L) [M_V]$')
        ax1.annotate(tt,[0.95,0.85],xycoords='axes fraction', ha='right', \
            size = 'medium' )
        ax1.set_xlim([-15,5])
        plt.gca().invert_xaxis()
        plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
        if k != 4:
            ax1.set_xticklabels(['' for i in range(10)])
            ax1.set_xlabel('')
        
        #ax2 = fig.add_subplot(1, 3, 2)
        ax2 = fig.add_subplot(gs2[k,1])
        print out.sfr_x[100]
        print out.sfr_x[-100]
        ax2.fill_between(ln2mv(out00.x), ln2mvpdf(table00[100]),\
                0, color=appa('#0000FF', alpha))
        ax2.fill_between(ln2mv(out00.x), ln2mvpdf(table00[-100]),\
                 0, color=appa('#00FF00',alpha))
        #print 'Make sure you know what SFRs these are'
        ax2.plot(ln2mv(out.x), ln2mvpdf(table[0]), color='blue')
        ax2.plot(ln2mv(out.x), ln2mvpdf(table[-1]), color='green')
        ax2.set_xlabel(r'$L_1 \, [M_V]$')
        ax2.set_ylabel(r'$p(L_1) \, [M_V]$')
        ax2.annotate(tt,[0.05,0.85],xycoords='axes fraction', ha='left', \
            size = 'medium' )
        ax2.set_xlim([-30, 0])
        plt.gca().invert_xaxis()
        plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
        if k != 4:
            ax2.set_xticklabels(['' for i in range(10)])
            ax2.set_xlabel('')
        
        #ax3 = fig.add_subplot(1, 3, 3)
        ax3 = fig.add_subplot(gs3[k,2])
        fill_between2(ax3, out.sfr_x, ln2mv(out.q95), ln2mv(out.q5),\
                       ln2mv(out00.q95), ln2mv(out00.q5),\
                       color=['#0000FF', '#FF0000'], alpha=[alpha,alpha])
        #ax3.fill_between(out.sfr_x, out.q5,out.q95, color='blue')
        #ax3.fill_between(out00.sfr_x, out00.q5,out00.q95, color=appa('#FF0000', alpha))
        ax3.set_xlabel(r'$\log \dot{M_{\star}} \ [M_{\odot} $ yr$^{-1}]$')
        ax3.set_ylabel(r'$L_1 [M_V]$')
        ax3.annotate(tt,[0.05,0.85],xycoords='axes fraction', ha='left', \
            size = 'medium' )
        ax3.set_ylim([-30,0])
        plt.gca().invert_yaxis()
    
        plt.gca().yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
        if k != 4:
            ax3.set_xticklabels(['' for i in range(10)])
            ax3.set_xlabel('')
        plt.tight_layout()
        print name
        #plt.show()
        k = k + 1
        return out, out00

    file00='dat/std.dat'
    d='dat/'

    k = 0
    fac = 0.75
    nrow = 5
    fig = plt.figure(figsize=(18*fac,4*fac*nrow))
    gs = gridspec.GridSpec(3, 1)
    gs1 = gridspec.GridSpec(nrow, 3)
    gs2 = gridspec.GridSpec(nrow, 3)
    gs3 = gridspec.GridSpec(nrow, 3)
    gs1.update(hspace=0)
    gs2.update(hspace=0)
    gs3.update(hspace=0)
    def_err = log(sqrt(2.512))
    #out, out00 = mk_plot(file00, d+'fc05.dat', r'$\mathcal{F}_c = 0.05$')
    out, out00 = mk_plot(file00, d+'fc001.dat', r'$\mathcal{F}_c = 0.001$', def_err)
    out, out00 = mk_plot(file00, d+'mmax12.dat', r'$m_{max}=10^{12} M_{\odot}$', def_err)
    out, out00 = mk_plot(file00, d+'mmax7.dat', r'$m_{max}=10^7 M_{\odot}$', def_err)
    out, out00 = mk_plot(file00, d+'mmax6.dat', r'$m_{max}=10^6 M_{\odot}$', def_err)
    out, out00 = mk_plot(file00, d+'mmin10000.dat', r'$m_{min}=10^4 M_{\odot}$', def_err)
    fig.savefig('grid1.eps')
    fig.show()

    k = 0
    fac = 0.75
    nrow = 5
    fig = plt.figure(figsize=(18*fac,4*fac*nrow))
    gs = gridspec.GridSpec(3, 1)
    gs1 = gridspec.GridSpec(nrow, 3)
    gs2 = gridspec.GridSpec(nrow, 3)
    gs3 = gridspec.GridSpec(nrow, 3)
    gs1.update(hspace=0)
    gs2.update(hspace=0)
    gs3.update(hspace=0)
    out, out00 = mk_plot(file00, d+'slope25.dat', r'$\alpha = -2.5$', def_err)
    out, out00 = mk_plot(file00, d+'slope15.dat', r'$\alpha = -1.5$', def_err)
    out, out00 = mk_plot(file00, d+'sfrerr.dat', r'$\log_{10} \sigma_{SFR} = 1.5$', def_err)
    out, out00 = mk_plot(file00, d+'obserr.dat', r'$\ln \sigma_L = 3$', 3.)
    out, out00 = mk_plot(file00, d+'dusty4.dat', r'$t_0 =0; t_1 = 5$', def_err)
    fig.savefig('grid2.eps')
    fig.show()
    out00, table00 = sfrl1_parse(file00)
    pdf00=out00.pdf_l

        
        

