
import math
import numpy as np
import matplotlib.pyplot as plt


phi, dphi = np.loadtxt('ens.dat').T

def meanvar(ith, Npeak, Ncosm):  #ith peak, num of peaks, num of cosmology pairs
	phi_dat  = [ phi[ith+i*Npeak] for i in range(Ncosm)]
	dphi_dat = [dphi[ith+i*Npeak] for i in range(Ncosm)]	
	dphi2    = [x*x for x in dphi_dat]

	phi_mean = sum(phi_dat)/Ncosm
	mean = sum(dphi_dat)/Ncosm
	var2 = sum(dphi2)/Ncosm - mean*mean
	return phi_mean, mean, math.sqrt(var2)

def f(x, nnu, nnu_fid):
	R  = 0.23*nnu/(1.+0.23*nnu)
	Rf = 0.23*nnu_fid/(1.+0.23*nnu_fid)
        R1 = 0.23/1.23
	R3 = 0.23*3/(1+0.23*3)
	A  = (R-Rf)/(R3-R1)
	if(x > trs):
		return A*(0.025+6.5e-4*x)
	else:
		return A*(0.025+6.5e-4*trs)/trs*x

Npeak = 8
Ncosm = 100

pphi = np.zeros(Npeak)
pdphi= np.zeros(Npeak)
errb = np.zeros(Npeak)

trs = 6*144*1.04e-2/math.pi
xfit = np.linspace(0, 8 , 100)
yfit = [-1*f(x, 1, 3) for x in xfit]

for ith in range(Npeak):
	pphi[ith], pdphi[ith], errb[ith] = meanvar(ith, Npeak, Ncosm)

alpha = np.linspace(2.,8.,100)
prob  = np.zeros(100)
for i in range(100):
	trs = alpha[i]*144*1.04e-2/math.pi
	model = [-1*f(x,1,3) for x in pphi]
	resid = [model[j]-pdphi[j] for j in range(len(model))]
	tmp   = [resid[j]**2/(2*errb[j]**2) for j in range(len(model))]
	prob[i] = sum(tmp)
	


plt.figure(0)
plt.errorbar(pphi, pdphi, yerr=errb, fmt= ' ', label = r'$\phi(N=1) -\phi(N=3)$')
plt.plot(xfit, yfit)
plt.plot((0,8),(0.042,0.042))
plt.legend(loc='best')
plt.xlabel(r'$k\times r_s(\eta(z=1000))/\pi$', fontsize = 20)
plt.ylabel(r'$\delta\phi/\pi$',fontsize=20)
plt.tight_layout()

plt.figure(1)
plt.plot(alpha, prob)
plt.show()

