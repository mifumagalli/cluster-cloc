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

sfr=[-2, -3, -1]
files=["sfrm2.txt", "sfrm3.txt","sfrm1.txt"]
i=0
#the same for all sfrs
fcs = 1e6*np.log(10)/(1e8-1e7)
#fcs2 = fcs*1611.81/qpow(100,1e9, -2., 0.6)
#multilpy by 6 looks good

#fcw.fall_cpp_wrapper(t0=0, t1=0.01, age_slope=-0.99,\
#                    obs_err=0.01, f_c=fcs, grid_out="slug_comp_rere.dat",\
#                   length=1e9,\
 #                  sfr_err=0.01, step=0.125/4, mmin=100, mmax=1e9)


#afcw.fall_cpp_wrapper(t0=0, t1=0.01, age_slope=-0.99,\
##                    obs_err=0.01, f_c=fcs2, grid_out="slug_comp.dat",\
#                   length=1e9, mmin=100,mmax=1e9\
#                    sfr_err=0.01, step=0.125/4)
fcw.fall_cpp_wrapper(t0=0, t1=0.01, age_slope=-0.99,\
                    obs_err=0.01, f_c=fcs, grid_out="slug_comp.dat",\
                   length=1e8,\
                   gamma_min = 2.438e19, sfr_err=0.01, step=0.125/4)

#fcs2 = 1e6*np.log(1000)/(1e9-1e8)
#fcw.fall_cpp_wrapper(t0=0, t1=0.01, age_slope=-0.99,\
#                    obs_err=0.01, f_c=fcs2, grid_out="slug_comp.dat",\
##                   length=1e9,\
 #                  sfr_err=0.01, step=0.125/4)

out,table = p.sfrl1_parse('slug_comp.dat')
#out2,table2 = p.sfrl1_parse('slug_comp2.dat')

#get names of tags
#out.dtype.names
sfrs = out.sfr_x
x = out.x
#get the sfr closest to \log sfr = -2
for i in xrange(3):
    index = np.argmin(abs(out.sfr_x -sfr[i]))
    this_pdf = table[index]
    #this_pdf2 = table2[index]
    print sfrs[index]
   
    plt.plot(p.ln2mv(x), p.ln2mvpdf(this_pdf), lw=3, color='blue', label='CLOC')
    #plt.plot(p.ln2mv(x), p.ln2mvpdf(this_pdf), lw=3, color='blue', label="Corrected")
    plt.xlim([-1,-18])
    f=open(files[i], 'r')
    slug = np.array(f.read().split("\n"))
    slug=slug[0:-1]
    slug = slug.astype(float)
    plt.hist(slug, normed=True, bins=40, histtype="stepfilled", color='grey',\
      edgecolor='none')
    plt.xlabel(r"$M_V$")
    plt.ylabel(r"$p(M_V)$")
    plt.title(r"$\log_{10}  \textrm{SFR} = "+repr(sfr[i])+"$")
    plt.legend(prop={'size':14}, frameon=False)
    plt.show()
    plt.savefig("slugcomp"+repr(i)+".eps")
    plt.clf()
    
    
