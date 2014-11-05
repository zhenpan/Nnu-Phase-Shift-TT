import template
import math
import numpy as np
import matplotlib.pyplot as plt


leq1, leq2, leq3 = np.loadtxt('leq.dat')

ell  = np.linspace(10,2500,100)
dell13 = [template.f(l, 3, 1, leq2) for l in ell]
dell35 = [template.f(l, 5, 3, leq2) for l in ell]

def bashinsky(nnu,nnu_fid):
	R  = 0.23*nnu/(1.+0.23*nnu)
	Rf = 0.23*nnu_fid/(1.+0.23*nnu_fid)
	dphi = (Rf-R)*0.1912*math.pi
	theta= 1.04136e-2
	return dphi/theta
print bashinsky(1,3), bashinsky(3,5)

plt.figure(0)
plt.plot(ell, dell13, 'r-', lw = 2, label=r'Template $\Delta\ell_{(N=1 - N=3)} = f(\ell)$')
plt.plot(ell, dell35, 'b--', lw = 2, label=r'Template $\Delta\ell_{(N=3-N=5)} = f(\ell)*0.57$')
plt.plot((10,2500),(12.7,12.7),'r-', lw = 2, label = r'Bashinsky $\Delta\ell_{(N=1-N=3)} = 12.7$')
plt.plot((10,2500),(7.3,7.3),'b--', lw = 2, label = r'Bashinsky $\Delta\ell_{(N=3-N=5)} = 7.3$')
plt.legend(loc='best')
plt.xlabel(r'$\ell$', fontsize = 20)
plt.ylabel(r'$\Delta\ell$', fontsize=20)
plt.savefig('fig2.pdf')
plt.show()
