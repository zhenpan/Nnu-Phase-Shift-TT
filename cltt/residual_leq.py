import template
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import camb4py

camb = camb4py.load('/home/zhenpan/camb/camb')

def run(As, ns, Tau, H0, omgb, omgc, nnu, Yp):
	result = camb(**{'get_scalar_cls':True, 'do_lensing':True, 'get_transfer':True, 'l_max_scalar': 6000,
		        'massless_neutrinos': nnu, 'massive_neutrinos':0, 'helium_fraction':Yp,
			'ombh2': omgb, 'omch2': omgc, 'hubble': H0, 're_optical_depth': Tau,
			'scalar_spectral_index(1)': ns,'scalar_amp(1)': As, 'pivot_scalar':0.05})
	Cls = result['lensed']
	Out = result['misc']
	return Cls, Out

def phase(phi, Delta):
	N_inte = len(phi)
	xnew = phi
	ynew = Delta
	raw  = [i for i in range(1, N_inte-1) if (ynew[i]-ynew[i-1])*(ynew[i]-ynew[i+1]) > 0. ]
	posi = [raw[i] for i in range(1,len(raw)) if(raw[i]-raw[i-1] > 50)]
	phi_m   = refine(phi, Delta, posi)
	Delta_m = [ynew[i] for i in posi]
	return  phi_m, Delta_m

def refine(phi, Delta, posi):
	phi_m = np.zeros([len(posi)])
	for i in range(len(posi)):
		j = posi[i]
		xnew = [phi[j-1], phi[j], phi[j+1]]
		ynew = [Delta[j-1], Delta[j], Delta[j+1]]
		p = np.polyfit(xnew, ynew, 2)
		phi_m[i] = -1*p[1]/(2.*p[0])
	return phi_m
	

def analyse(As, ns, Tau, H0, omgb, omgc, nnu, Yp):
	Cls, Out = run(As, ns, Tau, H0, omgb, omgc, nnu, Yp)
	k_D = float(Out['k_D(zstar) Mpc'])
	r_s = float(Out['r_s(zstar)/Mpc'])
	theta_s = float(Out['100*theta'])
	theta_D = float(Out['100*theta_D'])
	zeq = float(Out[' '])
	keq = math.sqrt(2*(omgb+omgc)*(1+zeq))/3000.
	leq = keq*r_s/(theta_s/100)
	con = 0.65
	l_D = con*math.pi/(theta_D/100.) #k_D*r_s/(theta_s/100.)
	ell = [l for l in range(2,4000)]
	tt  = [Cls[l-2,1]*math.exp(2*(l/l_D)**1.18) for l in ell]
	return ell, tt, theta_s, theta_D, l_D, leq




def rescale(k, Npeak): 

	ell1, tt1, ths1, thd1, ld1, leq1 = analyse(Amp1[k], Ns1[k], T_reio1[k], Hubble1[k], omgbh21[k], omgch21[k], Nnu1[k], Helium1[k])
	ell2, tt2, ths2, thd2, ld2, leq2 = analyse(Amp2[k], Ns2[k], T_reio2[k], Hubble2[k], omgbh22[k], omgch22[k], Nnu2[k], Helium2[k])
	
	nnu1 = Nnu1[k]
	nnu2 = Nnu2[k]

	ell1 = [l + template.f(l, nnu1, nnu2, leq2) for l in ell1]

	ell1_m, tt1_m = phase(ell1, tt1)
	ell2_m, tt2_m = phase(ell2, tt2)

	plt.figure(0)
	plt.plot(ell1, tt1, label = r'${N_{\rm eff} = 1}$')
	plt.plot(ell2, tt2, label = r'${N_{\rm eff} = 3}$')
	plt.plot(ell1_m, tt1_m, 'o', ell2_m, tt2_m, 'o')
	plt.legend(loc='lower right')
	plt.xlabel(r'$\ell$',fontsize=20)
	plt.ylabel(r'$C^{TT}_{\ell\ell}$',fontsize=20)
	plt.tight_layout()
#	plt.show()

	ell_m = [ell1_m[i] for i in range(Npeak)]
	dell  = [ell1_m[i]-ell2_m[i] for i in range(Npeak)]
	print k, nnu1, nnu2, leq1, leq2

	plt.figure(1)
	plt.plot(ell_m, dell)
	plt.plot(ell_m, dell, 'o', label = r'$\phi(N=1) -\phi(N_{\rm fid})$')
	plt.plot((0, 3000),(0.,0.), '-')
	plt.xlim(0,3000)
	plt.legend(loc='best')
	plt.xlabel(r'$\ell$', fontsize = 20)
	plt.ylabel(r'$\delta\ell$',fontsize=20)
	plt.title(r'$\ell$')
	plt.tight_layout()
#	plt.show()

	tt1  = [tt1[i]*template.g(ell1[i], nnu1, nnu2,leq2) for i in range(len(ell1))]


	plt.figure(2)
	plt.plot(ell1, tt1, label = r'${N_{\rm eff} = 1}$')
	plt.plot(ell2, tt2, label = r'${N_{\rm eff} = 3}$')
	plt.legend(loc='lower right')
	plt.xlabel(r'$\ell$',fontsize=20)
	plt.legend(loc='lower right')
	plt.xlabel(r'$\ell$',fontsize=20)
	plt.ylabel(r'$C^{TT}_{\ell\ell}$',fontsize=20)
	plt.tight_layout()
	plt.show()

	return dell


Amp1, Ns1 , T_reio1, Hubble1, omgbh21, omgch21, Nnu1, Helium1 = np.loadtxt('params1a').T
Amp2, Ns2 , T_reio2, Hubble2, omgbh22, omgch22, Nnu2, Helium2 = np.loadtxt('params1b').T
Npeak = 18
Nsamp = 100
res = np.zeros((Nsamp, Npeak))

for k in range(Nsamp):
	error = rescale(k, Npeak) 
	for m in range(Npeak):
		res[k][m] = error[m]

np.savetxt('res1.dat', res)	
