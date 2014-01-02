import fall_cpp_wrapper as fcw
import numpy as np
import parse as p
import matplotlib.pyplot as plt
from matplotlib import rc
def qpow(min1, max1, slope, q):
    return (q*(max1**(slope+1) - min1**(slope+1)) + min1**(slope+1))**(1./(slope+1))

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


fcw.fall_cpp_wrapper(t0=0, t1=1.0, age_slope=-0.9,\
                   obs_err=0.005, f_c=0.01, grid_out="slug_comp.dat",\
                   length=1e8,mmin=500, mmax=1e9,\
                   gamma_min = 2.438e19, sfr_err=0.005, step=0.125/4)

out,table = p.sfrl1_parse('slug_comp.dat')

mcdata = open("/Users/rdasilva/Dropbox/sfr_l1_plots/montecarlo/mctest.txt",'r')
mcdata = mcdata.read().split("\n")[0:-1]
mcdata = np.array(mcdata).astype(float)
mcdata=mcdata[mcdata<-3]

#get names of tags
#out.dtype.names
sfrs = out.sfr_x
x = out.x
#get the sfr closest to \log sfr = -1
index = np.argmin(abs(out.sfr_x -(-1)))
this_pdf = table[index]
   
plt.xlim([-1,-18])

plt.hist(mcdata, normed=True, bins=80*4, histtype="stepfilled", color='grey',\
  edgecolor='none')

plt.plot(p.ln2mv(x), p.ln2mvpdf(this_pdf), lw=3, color='blue')
plt.xlabel(r"$M_V$")
plt.ylabel(r"$p(M_V)$")
plt.show()
plt.savefig("compmc.eps")
plt.clf()


