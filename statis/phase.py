import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def phase(kmode, phi, Delta, kD):
	N_inte = 40000
	kmode_c= [kmode[i] for i in range(1,phi.size) if (phi[i]-phi[i-1]) > 1e-4]
	phi_c  = [phi[i]   for i in range(1,phi.size) if (phi[i]-phi[i-1]) > 1e-4]
	Delta_c= [Delta[i] for i in range(1,phi.size) if (phi[i]-phi[i-1]) > 1e-4]
	Delta_c= [Delta[i]*math.exp((kmode[i]/kD)**2) for i in range(1,phi.size) if (phi[i]-phi[i-1]) > 1e-4]
	f    = interp1d(phi_c, Delta_c, kind='cubic')
	xnew = np.linspace(phi_c[0], phi_c[-1], N_inte)
	ynew = f(xnew)
	posi    = [i for i in range(1, N_inte-1) if (ynew[i]-ynew[i-1])*(ynew[i]-ynew[i+1]) > 0. ]
	phi_m   = [xnew[i] for i in posi]
	Delta_m = [ynew[i] for i in posi]	
	return  phi_c, Delta_c, phi_m, Delta_m

def f(x, nnu, nnu_fid):
	R  = 0.23*nnu/(1.+0.23*nnu)
	Rf = 0.23*nnu_fid/(1.+0.23*nnu_fid)
        R1 = 0.23/1.23
	R3 = 0.23*3/(1+0.23*3)
	A  = (R-Rf)/(R3-R1)
	trs= 2.8
	if(x > trs):
		return A*(0.025+6.5e-4*x)
	else:
		return A*(0.025+6.5e-4*trs)/trs*x
#	return A*(0.028/2.84*x/(math.exp(10*(x-2.84))+1) + (0.026+6.5e-4*x)/(math.exp(10*(2.84-x))+1.))

def scale(nnu, nnu_fid):
	R  = 0.23*nnu/(1.+0.23*nnu)
	Rf = 0.23*nnu_fid/(1.+0.23*nnu_fid)
        R1 = 0.23/1.23
	R3 = 0.23*3/(1+0.23*3)
	return  (R3-R1)/(Rf-R)

def process(ith, kmode, phi, Delta, k_D):
	start = ith*200
	end   = (ith+1)*200
	phi_c, Delta_c, phi_m, Delta_m = phase(kmode[start:end], phi[start:end], Delta[start:end], k_D[ith])
	return phi_c, Delta_c, phi_m, Delta_m


def stat(ith):
	phi1_c, Delta1_c, phi1_m, Delta1_m = process(ith, kmode1, phi1, Delta1, k_D1)
	phi2_c, Delta2_c, phi2_m, Delta2_m = process(ith, kmode2, phi2, Delta2, k_D2)

	norm = scale(Nnu1[ith], Nnu2[ith])
	dphi = [norm*(phi1_m[i]-phi2_m[i]) for i in range(len(phi2_m))]

#	plt.figure(1)
#	plt.plot(phi1_c, Delta1_c)
#	plt.plot(phi2_c, Delta2_c)

	return phi1_m, dphi


Nnu1 = np.loadtxt('params1a', usecols=(6,)).T
Nnu2 = np.loadtxt('params1b', usecols=(6,)).T

k_D1, k_D2     = np.loadtxt('difu.dat').T
kmode1, phi1, Delta1 = np.loadtxt('data1').T
kmode2, phi2, Delta2 = np.loadtxt('data2').T

xfit = np.linspace(0, 8 , 100)
yfit = [-1*f(x, 1, 3) for x in xfit]

phi  = []
dphi = []

for i in range(100):
	print i, Nnu1[i], Nnu2[i]
	phi_dat, dphi_dat = stat(i)
	phi  = phi + phi_dat
	dphi = dphi + dphi_dat 
out = np.column_stack((phi, dphi))
np.savetxt('ens.dat', out )

plt.figure(0)
plt.plot(phi, dphi, 'o', label = r'$\phi(N=1) -\phi(N=3)$')
plt.plot(xfit, yfit)
plt.plot((0,8),(0.042,0.042))
plt.legend(loc='best')
plt.xlabel(r'$k\times r_s(\eta(z=1000))/\pi$', fontsize = 20)
plt.ylabel(r'$\delta\phi/\pi$',fontsize=20)
plt.tight_layout()
plt.savefig('fig2.pdf')
plt.show()

