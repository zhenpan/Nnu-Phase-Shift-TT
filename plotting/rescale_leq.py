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
	ell = [l for l in range(2,3000)]
	tt  = [Cls[l-2,1]*math.exp(2*(l/l_D)**1.18) for l in ell]
	return ell, tt, theta_s, theta_D, l_D, leq



As, ns , Tau, H0, omgb, omgc, nnu1, Yp = np.loadtxt('params', usecols= (1,))
ell, tt1, ths1, thd1, ld1, leq1 = analyse(As, ns, Tau, H0, omgb, omgc, nnu1, Yp)

As, ns , Tau, H0, omgb, omgc, nnu2, Yp = np.loadtxt('params', usecols= (2,))
ell, tt2, ths2, thd2, ld2, leq2 = analyse(As, ns, Tau, H0, omgb, omgc, nnu2, Yp)

As, ns , Tau, H0, omgb, omgc, nnu3, Yp = np.loadtxt('params', usecols= (3,))
ell, tt3, ths3, thd3, ld3, leq3 = analyse(As, ns, Tau, H0, omgb, omgc, nnu3, Yp)
	

nnu_fid = 3
leq_fid = leq2

ell1 = [l + template.f(l, nnu1, nnu_fid, leq_fid) for l in ell]
ell2 = [l + template.f(l, nnu2, nnu_fid, leq_fid) for l in ell]
ell3 = [l + template.f(l, nnu3, nnu_fid, leq_fid) for l in ell]

ell1_m, tt1_m = phase(ell1, tt1)
ell2_m, tt2_m = phase(ell2, tt2)
ell3_m, tt3_m = phase(ell3, tt3)

plt.figure(0)
plt.plot(ell1, tt1, label = r'${N_{\rm eff} = 1}$')
plt.plot(ell2, tt2, label = r'${N_{\rm eff} = 3}$')
plt.plot(ell3, tt3, label = r'${N_{\rm eff} = 5}$')
plt.plot(ell1_m, tt1_m, 'o', ell2_m, tt2_m, 'o', ell3_m, tt3_m, 'o')
plt.legend(loc='lower right')
plt.xlabel(r'$\ell$',fontsize=20)
plt.ylabel(r'$C^{TT}_{\ell\ell}$',fontsize=20)
plt.tight_layout()
plt.savefig('fig2.pdf')
plt.show()

Npeak  = min(len(ell1_m),len(ell2_m),len(ell3_m))
dell12 = [ell1_m[i]-ell2_m[i] for i in range(Npeak)]
dell23 = [ell2_m[i]-ell3_m[i] for i in range(Npeak)]
ell_m  = [ell1_m[i] for i in range(Npeak)]
dell   = [ell1_m[i]-ell2_m[i] for i in range(Npeak)]

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
plt.savefig('fig3.pdf')
plt.show()

tt1  = [tt1[i]*g(ell[i], nnu1, nnu_fid, leq_fid) for i in range(len(ell))]
tt2  = [tt2[i]*g(ell[i], nnu2, nnu_fid, leq_fid) for i in range(len(ell))]
tt3  = [tt3[i]*g(ell[i], nnu3, nnu_fid, leq_fid) for i in range(len(ell))]

tt1_m  = [tt1_m[i]*g(ell1_m[i], nnu1, nnu_fid, leq_fid) for i in range(len(ell3_m))]
tt2_m  = [tt2_m[i]*g(ell2_m[i], nnu2, nnu_fid, leq_fid) for i in range(len(ell3_m))]
tt3_m  = [tt3_m[i]*g(ell3_m[i], nnu3, nnu_fid, leq_fid) for i in range(len(ell3_m))]


def interp(ell, tt, samp):
	tts  = np.zeros([len(samp)])
	j = 0
	for i in range(len(samp)):
		while(ell[j] < samp[i]):
			j = j+1
		xnew = [ell[j-1], ell[j], ell[j+1]]
		ynew = [tt[j-1], tt[j], tt[j+1]]
		z    = np.polyfit(xnew, ynew, 2)
		p    = np.poly1d(z)
		tts[i]= p(samp[i])
	return tts

samp = np.linspace(1, 2500, num =40) 
tts1 = interp(ell1, tt1, samp)
tts2 = interp(ell2, tt2, samp)
tts3 = interp(ell3, tt3, samp)

plt.figure(2)
plt.plot(ell1, tt1, label = r'${N_{\rm eff} = 1}$')
plt.plot(ell2, tt2, label = r'${N_{\rm eff} = 3}$')
plt.plot(ell3, tt3, label = r'${N_{\rm eff} = 5}$')
plt.legend(loc='lower right')
plt.xlabel(r'$\ell$',fontsize=20)
plt.legend(loc='lower right')
plt.xlabel(r'$\ell$',fontsize=20)
plt.ylabel(r'$C^{TT}_{\ell\ell}$',fontsize=20)
plt.tight_layout()
plt.savefig('fig3.pdf')
plt.show()

damp12 = [tts1[i]/tts2[i] for i in range(len(samp))]
damp23 = [tts2[i]/tts3[i] for i in range(len(samp))]

plt.figure(3)
plt.plot(samp, damp12, label = r'${\rm Amp}_{N=1}/{\rm Amp}_{N=fid}$')
plt.plot(samp, damp23, label = r'${\rm Amp}_{N=fid}/{\rm Amp}_{N=5}$')
plt.plot(samp, damp12,'o', samp, damp23, 'o')
plt.legend(loc='best')
plt.xlabel(r'$\ell$',fontsize=20)
plt.tight_layout()
plt.show()

