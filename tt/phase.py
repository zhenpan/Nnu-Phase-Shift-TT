
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
	Cls = result['scalar']
	Out = result['misc']
	return Cls, Out


def phase(phi, Delta):
	N_inte = len(phi)
	xnew = phi
	ynew = Delta
	raw  = [i for i in range(1, N_inte-1) if (ynew[i]-ynew[i-1])*(ynew[i]-ynew[i+1]) > 0. and (ynew[i]-ynew[i-5])*(ynew[i]-ynew[i+5]) > 0.]
	posi = [raw[i] for i in range(1,len(raw)) if(raw[i]-raw[i-1] > 10)]
	phi_m, Delta_m  = refine(phi, Delta, posi)
	return  phi_m, Delta_m

def refine(phi, Delta, posi):
	phi_m = np.zeros([len(posi)])
	Delta_m = np.zeros([len(posi)])
	for i in range(len(posi)):
		j = posi[i]
		xnew = [phi[j-1], phi[j], phi[j+1]]
		ynew = [Delta[j-1], Delta[j], Delta[j+1]]
		z    = np.polyfit(xnew, ynew, 2)
		p    = np.poly1d(z)
		phi_m[i]   = -1*z[1]/(2.*z[0])
		Delta_m[i] = p(phi_m[i]) 
	return phi_m, Delta_m 
	
	

def analyse(As, ns, Tau, H0, omgb, omgc, nnu, Yp):
	Cls, Out = run(As, ns, Tau, H0, omgb, omgc, nnu, Yp)
	k_D = float(Out['k_D(zstar) Mpc'])
	r_s = float(Out['r_s(zstar)/Mpc'])
	theta_s = float(Out['100*theta'])
	theta_D = float(Out['100*theta_D'])
	zeq = float(Out[' '])
	con = 0.65
	l_D = con*math.pi/(theta_D/100.)              #k_D*r_s/(theta_s/100.)
	keq = math.sqrt(2*(omgb+omgc)*(1+zeq))/3000.
	leq = keq*r_s/(theta_s/100)
	norm= 1. #/Cls[700,1]
	ell = [l for l in range(2,4000)]
	tt  = [norm*Cls[l-2,1]*math.exp(2*(l/l_D)**1.18) for l in ell]
	return ell, tt, theta_s, theta_D, zeq, l_D, leq


As, ns , Tau, H0, omgb, omgc, nnu, Yp = np.loadtxt('params', usecols= (1,))
ell, tt1, ths1, thd1, zeq1, ld1, leq1 = analyse(As, ns, Tau, H0, omgb, omgc, nnu, Yp)

As, ns , Tau, H0, omgb, omgc, nnu, Yp = np.loadtxt('params', usecols= (2,))
ell, tt2, ths2, thd2, zeq2, ld2, leq2 = analyse(As, ns, Tau, H0, omgb, omgc, nnu, Yp)

As, ns , Tau, H0, omgb, omgc, nnu, Yp = np.loadtxt('params', usecols= (3,))
ell, tt3, ths3, thd3, zeq3, ld3, leq3 = analyse(As, ns, Tau, H0, omgb, omgc, nnu, Yp)
	
print  '100*theta_star:', ths1, ths2, ths3
print  '100*theta_D:', thd1, thd2, thd3
print  'z_eq:', zeq1, zeq2, zeq3
print  'l_D:', ld1,  ld2,  ld3
print 'l_eq:', leq1, leq2, leq3

ell1_m, tt1_m = phase(ell, tt1)
ell2_m, tt2_m = phase(ell, tt2)
ell3_m, tt3_m = phase(ell, tt3)

#X = [i for i in range(len(ell1_m))]
#p1 = np.polyfit(X, ell1_m,1)
#p2 = np.polyfit(X, ell2_m,1)
#p3 = np.polyfit(X, ell3_m,1)
#print p1, p2, p3
#Y1 = [p1[0]*x+p1[1]-ell1_m[x] for x in X]
#Y2 = [p2[0]*x+p2[1]-ell2_m[x] for x in X]
#Y3 = [p3[0]*x+p3[1]-ell3_m[x] for x in X]
#plt.plot(X, Y1, 'o', X, Y2, 'o', X, Y3, 'o')
#plt.plot(X, Y1,  X, Y2,  X, Y3 )
#plt.show()

plt.figure(0)
plt.plot(ell, tt1, label = r'${N_{\rm eff} = 1}$')
plt.plot(ell, tt2, label = r'${N_{\rm eff} = 3}$')
plt.plot(ell, tt3, label = r'${N_{\rm eff} = 5}$')
plt.plot(ell1_m, tt1_m, 'o')
plt.plot(ell2_m, tt2_m, 'o')
plt.plot(ell3_m, tt3_m, 'o')
plt.legend(loc='lower right')
plt.xlabel(r'$\ell$',fontsize=20)
plt.ylabel(r'$C^{TT}_{\ell\ell}$',fontsize=20)
plt.tight_layout()
plt.savefig('fig2.pdf')
plt.show()
	
dell12 = [ell1_m[i]-ell2_m[i] for i in range(len(ell1_m))]
dell23 = [ell2_m[i]-ell3_m[i] for i in range(len(ell1_m))]
out = np.column_stack((ell1_m, dell12, dell23))
np.savetxt('dell.dat', out)

elleq = [leq1, leq2, leq3]
np.savetxt('leq.dat', elleq)

plt.figure(1)
plt.plot(ell1_m, dell12 )
plt.plot(ell1_m, dell23 )
plt.plot(ell1_m, dell12, 'o', label = r'$\phi(N=1) -\phi(N=3)$')
plt.plot(ell1_m, dell23, 'o', label = r'$\phi(N=3) -\phi(N=5)$')
plt.legend(loc='upper left')
plt.xlabel(r'$\ell$', fontsize = 20)
plt.ylabel(r'$\delta\ell$',fontsize=20)
plt.title(r'$\ell$')
plt.tight_layout()
plt.savefig('fig3.pdf')
plt.show()

