
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


def analyse(As, ns, Tau, H0, omgb, omgc, nnu, Yp):
	Cls, Out = run(As, ns, Tau, H0, omgb, omgc, nnu, Yp)
	k_D = float(Out['k_D(zstar) Mpc'])
	r_s = float(Out['r_s(zstar)/Mpc'])
	theta_s = float(Out['100*theta'])
	theta_D = float(Out['100*theta_D'])
	con = 0.65
	l_D = con*math.pi/(theta_D/100.) #k_D*r_s/(theta_s/100.)
	print l_D
	ell = [l for l in range(2,3000)]
	tt  = [Cls[l-2,1]*math.exp(2*(l/1267.1)**1.18) for l in ell]  #lD_fid = 1267.1
	return ell, tt, theta_s, theta_D, l_D

def f(x, nnu, nnu_fid): 
	R  = 0.23*nnu/(1.+0.23*nnu)
	Rf = 0.23*nnu_fid/(1.+0.23*nnu_fid)
        R1 = 0.23/1.23
	R3 = 0.23*3/(1+0.23*3)
	A  = (R-Rf)/(R3-R1)
	return A*(8.22 + 8.75e-4*x-4.2*math.exp(-0.1*(x-600.))/(math.exp(-0.1*(x-600.))+1)-4/(math.exp(0.05*(x-200))+1))


As, ns , Tau, H0, omgb, omgc, nnu1, Yp = np.loadtxt('params', usecols= (1,))
ell1, tt1, ths1, thd1, ld1 = analyse(As, ns, Tau, H0, omgb, omgc, nnu1, Yp)

As, ns , Tau, H0, omgb, omgc, nnu2, Yp = np.loadtxt('params', usecols= (2,))
ell2, tt2, ths2, thd2, ld2 = analyse(As, ns, Tau, H0, omgb, omgc, nnu2, Yp)

As, ns , Tau, H0, omgb, omgc, nnu3, Yp = np.loadtxt('params', usecols= (3,))
ell3, tt3, ths3, thd3, ld3 = analyse(As, ns, Tau, H0, omgb, omgc, nnu3, Yp)
	
nnu_fid = 3
nell  = 200

tt1  = [tt1[i]/1000. for i in range(len(ell1))] #amp norm
tt2  = [tt2[i]/1000. for i in range(len(ell2))]
tt3  = [tt3[i]/1000. for i in range(len(ell3))]

norm1 = tt2[nell]/tt1[nell]
norm2 = 1.
norm3 = tt2[nell]/tt3[nell]
tt1n = [tt1[i]*norm1 for i in range(len(ell1))]
tt2n = [tt2[i]*norm2 for i in range(len(ell2))]
tt3n = [tt3[i]*norm3 for i in range(len(ell3))]

#ell1c = [l + f(l, nnu1, nnu_fid) for l in ell1] #phase shift correction
#ell2c = [l + f(l, nnu2, nnu_fid) for l in ell2]
#ell3c = [l + f(l, nnu3, nnu_fid) for l in ell3]


As, ns , Tau, H0, omgb, omgc, nnu1, Yp = np.loadtxt('params_bbn', usecols= (1,))
ell_bbn1, tt_bbn1, ths1, thd1, ld1 = analyse(As, ns, Tau, H0, omgb, omgc, nnu1, Yp)

As, ns , Tau, H0, omgb, omgc, nnu1, Yp = np.loadtxt('params_bbn', usecols= (2,))
ell_bbn2, tt_bbn2, ths1, thd1, ld1 = analyse(As, ns, Tau, H0, omgb, omgc, nnu1, Yp)

As, ns , Tau, H0, omgb, omgc, nnu1, Yp = np.loadtxt('params_bbn', usecols= (3,))
ell_bbn3, tt_bbn3, ths1, thd1, ld1 = analyse(As, ns, Tau, H0, omgb, omgc, nnu1, Yp)

tt_bbn1  = [tt_bbn1[i]/1000. for i in range(len(ell_bbn1))] #amp norm
tt_bbn2  = [tt_bbn2[i]/1000. for i in range(len(ell_bbn2))]
tt_bbn3  = [tt_bbn3[i]/1000. for i in range(len(ell_bbn3))]


X, Y, Yerr = np.loadtxt('planck.dat', delimiter=',', usecols=(0,3,4)).T
Y = [math.exp(2*(X[i]/1267.1)**1.18)/1000*Y[i] for i in range(X.size)]
Yerr = [math.exp(2*(X[i]/1267.1)**1.18)/1000*Yerr[i] for i in range(X.size)]


plt.figure(0)
f, (c,a, b) = plt.subplots(3,sharex=True, sharey=True, figsize=(12,6))
f.set_size_inches(3,1.5)
c.plot(ell_bbn1, tt_bbn1, label = r'${N_{\rm eff} = 1}$')
c.plot(ell_bbn2, tt_bbn2, label = r'${N_{\rm eff} = 3}$')
c.plot(ell_bbn3, tt_bbn3, label = r'${N_{\rm eff} = 5}$')
c.errorbar(X,Y, yerr = Yerr, fmt = ' ', color = 'k')
c.legend(loc='lower right', prop={'size':12})
a.plot(ell1, tt1, label = r'${N_{\rm eff} = 1}$')
a.plot(ell2, tt2, label = r'${N_{\rm eff} = 3}$')
a.plot(ell3, tt3, label = r'${N_{\rm eff} = 5}$')
a.errorbar(X,Y, yerr = Yerr, fmt = ' ', color = 'k')
b.plot(ell1, tt1n)
b.plot(ell2, tt2n)
b.plot(ell3, tt3n)
b.errorbar(X,Y, yerr = Yerr, fmt = ' ', color = 'k', label = r'Planck Data')
b.legend(loc='lower right', prop={'size':12})
f.subplots_adjust(hspace=0)
plt.xlim(0,2500)
plt.xlabel(r'$\ell$',fontsize=20)
f.text(0.05,0.5, r'$\exp(2*(\ell/\ell_D^{\rm fid})^{1.18})\times\ell(\ell+1)C^{TT}_{\ell}/(2\pi) \ [10^3 \mu k^2]$',fontsize=15, ha='center', va='center', rotation='vertical')
f.text(0.5, 0.7, r'$\omega_b, \ z_{\rm eq}, \ \theta_s \ {\rm fixed}$ ', ha='center', va='center')
f.text(0.5, 0.44, r'$\omega_b, \ z_{\rm eq}, \ \theta_s, \ \theta_D \ {\rm fixed}$ ', ha='center', va='center')
f.text(0.5, 0.18,  r'$C_{200}(N_{\rm eff}) = C_{200}^{\rm best}$', ha='center', va='center')
f.text(0.5, 0.14, r'$\omega_b, \ z_{\rm eq}, \ \theta_s, \ \theta_D \ {\rm fixed}$ ', ha='center', va='center')
plt.savefig('fig1.pdf')
plt.show()

