#a#
import math
import numpy as np
from scipy.integrate import quad



def rs(x,omgb, omgc,nnu, H0):
	omgg  = 2.474e-5    	       
	omgr  = (2.474 + 0.562*nnu)*1.e-5    	       
	omgm  = omgb + omgc
	omgl  = (H0/100.)**2 - omgm - omgr
	ratio = 0.75*(omgb/omgg)*x	
	return 3000./math.sqrt(omgr + omgm*x + omgl*x**4)/math.sqrt(3*(1+ratio))


def thermo(x, omgb, omgc, H0,nnu, Yp):
	omgg  = 2.474e-5    	       
	omgr  = (2.474 + 0.562*nnu)*1.e-5    	       
	omgm  = omgb + omgc
	omgl  = (H0/100.)**2 - omgm - omgr
	ratio = 0.75*(omgb/omgg)*x	
	const = (ratio*ratio+16./15.*(1+ratio))/6/(1+ratio)**2
	nsiga3 = 2.3e-5*omgb*(1-Yp/2.)
	return 3000/math.sqrt(omgr + omgm*x + omgl*x**4)*(x*x)*(const/nsiga3)

def k_difu(omgb, omgc, H0,nnu, Yp):  # kD of pre-recombination
	r_d2 = quad(thermo, 0., a_rec, args=(omgb, omgc, H0, nnu, Yp))[0]
	return 1./math.sqrt(r_d2)

Amp, Ns, T_reio, Hubble1, omgbh21, omgch21, Nnu1, Helium1 = np.loadtxt('params1').T
Amp, Ns, T_reio, Hubble2, omgbh22, omgch22, Nnu2, Helium2 = np.loadtxt('params2').T

a_rec = 1./1000
k_D1 = np.zeros(len(Amp))
k_D2 = np.zeros(len(Amp))

for m in range(len(Amp)):
	k_D1[m]  = k_difu(omgbh21[m], omgch21[m], Hubble1[m], Nnu1[m], Helium1[m])
	k_D2[m]  = k_difu(omgbh22[m], omgch22[m], Hubble2[m], Nnu2[m], Helium2[m])

out = np.column_stack((k_D1, k_D2))
np.savetxt('difu.dat', out)
