#! /usr/bin/env python
import h5py
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, AutoMinorLocator,MultipleLocator

plt.rcParams['text.latex.preamble'] = r'\usepackage{mathrsfs}'
plt.rcParams['xtick.minor.visible'] = r'True'
plt.rcParams['ytick.direction'] = r'in'
plt.rcParams['xtick.direction'] = r'in'

plt.rcParams['ytick.right'] = r'True'
plt.rcParams['xtick.top'] = r'True'
plt.rcParams.update({"text.usetex": True})

nbin=24
cgnc='#9467BD'
cgns='#1F77B4'
#2CA02C
cgnc1="#2CA02C"
cana1='#FF7F0E'
cgnc_lc="#6798BD"
colorstar=cgnc
colorsbh=cgnc1
colorns=cgns
mzstar=5
mzsbh=5
#mzbbh=5
mzns=4
stylestar='.-'
stylesbh='o'
#stylebbh='-o'
stylens='^'
stylestar_ana='-'
stylesbh_ana='--'
#stylebbh='-o'
stylens_ana=':'
rmin=np.log10(3.1*206264/1e5/2.)
rb=np.log10(3.1*206264*2*2)
rh=np.log10(3.1*206264)
nh=np.log10(2e4/206264.**3)
lpc3=np.log10(206264)*3


def set_xyaxis(ax=None, xmajorstep=None,xminorstep=None, \
			   ymajorstep=None,yminorstep=None, \
			   xmajorticks=5,ymajorticks=5,xminorticks=5,\
			   yminorticks=5,scifmt=False, xvisible=True, yvisible=True):
	plt.minorticks_on()
	if(ax==None):
		ax=plt.gca()
	if(scifmt==True):
		set_sci_format(ax,'both')
	if(xmajorstep == None):
		ax.xaxis.set_major_locator(MaxNLocator(nbins=xmajorticks))
	else:
		ax.xaxis.set_major_locator(MultipleLocator(base=xmajorstep))
	if(xminorstep == None):
		ax.xaxis.set_minor_locator(AutoMinorLocator(n=xminorticks))
	else:
		ax.xaxis.set_minor_locator(MultipleLocator(base=xminorstep))
		
	if(ymajorstep == None):	
		ax.yaxis.set_major_locator(MaxNLocator(nbins=ymajorticks))
	else:
		ax.yaxis.set_major_locator(MultipleLocator(base=ymajorstep))		
	if(yminorstep == None):
		ax.yaxis.set_minor_locator(AutoMinorLocator(n=yminorticks))
	else:	
		ax.yaxis.set_minor_locator(MultipleLocator(base=yminorstep))
	if(yvisible==False):
		ax.yaxis.set_ticklabels([])
	if(xvisible==False):
		ax.xaxis.set_ticklabels([])

def get_scatters_main(fdir, sm,sn,sb,cr=False):
	nums=1; nume=5; ncol=nbin
	x=np.zeros([nume-nums, ncol])
	y=np.zeros([nume-nums, ncol])
	for j in range(nums, nume):
		f=h5py.File(fdir+str(sm)+"_"+str(sn-j)+".hdf5", 'r' )
		#print(sn,j-nums)
		x[j-nums,:]=f[sb]['   X']
		y[j-nums,:]=f[sb]['  FX']
	if(cr):
		mean=np.mean(np.log10(y)+lpc3,axis=0)
		err=np.var(np.log10(y)+lpc3,axis=0)**0.5
		r=x[0,:]-rh
	else:
		mean=np.mean(np.log10(y),axis=0)
		err=np.var(np.log10(y),axis=0)**0.5
		r=x[0,:]
	return r,mean,err

def plot_one(ax,fdir, sb ,sn, sp, style, label):
	f=h5py.File(fdir+str(sn)+"_"+str(sp)+".hdf5", 'r' )
	x=f[sb]['   X']
	y=f[sb]['  FX']
	ax.plot(x,y,style,label=label)
	return x,y

	

plt.figure(figsize=(8,4))
plt.clf()

sx=1;sy=2

ax=plt.subplot(sx,sy,1)
for i in [1,2,3,4,5,10]:
	r, mean, err=get_scatters_main("../output/ecev/dms/dms_", i,10,  '1/star/fgx/')	
	ax.errorbar(r, mean, yerr=err,fmt=stylestar,label=f'{i/10.:.1f}'+' $T_{\\rm rlx}$',mfc='w', markersize=mzstar)

ax.legend(loc='upper left',ncol=2)

#ax.legend(loc="lower right", ncol=1)
ax.set_xlim(np.log10(0.05),5.2)
ax.set_ylim(-0.2,1.2)
#ax.set_yscale("log")
ax.set_xlabel("log $x$",fontsize=18)
ax.set_ylabel("log $\\bar g(x)$",fontsize=18)
set_xyaxis(xmajorstep=1,ymajorstep=1)

#######################################################################################################################
ax=plt.subplot(sx,sy,2)
#plot_scatters(ax, "../../../version_1.7_1D/model/nolc_1comp/output/ecev/dms/dms_burn_in_", 50,"1/star/fden",\
#	stylestar,cana1,mz=mzstar,mfc='w')
#plot_scatters(ax, "../../../version_1.7_1D/model/nolc_1comp/output/ecev/dms/dms_burn_in_", 50,"1/star/fden",\
#	"o-",cgnc1,mz=mzstar)	
for i in [1,2,3,4,5,10]:
	r, mean, err=get_scatters_main("../output/ecev/dms/dms_", i,10,  '1/star/fden/',cr=True)	
	ax.errorbar(r, mean, yerr=err,fmt=stylestar,label=str(i),mfc='w', markersize=mzstar)	

#ax.set_xscale("log")
ax.set_xlim(rmin-rh,0)
ax.set_ylim(3,13.5)
#ax.set_yscale("log")
ax.set_xlabel("log $r/r_h$",fontsize=18)
ax.set_ylabel("log $n(r)$ (pc$^{-3})$",fontsize=18)
set_xyaxis(xmajorstep=1.0,ymajorstep=2.0)

plt.tight_layout()
plt.savefig("fig_gE_1_lc.pdf")
