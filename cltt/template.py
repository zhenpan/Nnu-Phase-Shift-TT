import math


def f(x, nnu, nnu_fid, leq_fid):
	R  = 0.23*nnu/(1.+0.23*nnu)
	Rf = 0.23*nnu_fid/(1.+0.23*nnu_fid)
        R1 = 0.23/1.23
	R3 = 0.23*3/(1+0.23*3)
	A  = (R-Rf)/(R3-R1)
	trs = 5.5*leq_fid
	if(x > trs):
		return A*(8.43+6.5e-4*x)
	else: 
		return A*(8.43+6.5e-4*trs)/trs*x
	#return A*((8.43+6.5e-4*x)/(math.exp(-0.05*(x-trs))+1) + 1.47e-2*x/(math.exp(0.05*(x-trs))+1))


def g(x, nnu, nnu_fid, leq_fid):
	R  = 0.23*nnu/(1.+0.23*nnu)
	Rf = 0.23*nnu_fid/(1.+0.23*nnu_fid)
        R1 = 0.23/1.23
	R3 = 0.23*3/(1+0.23*3)
	B  = (R-Rf)/(R3-R1)
	trs= 2*leq_fid 
	tmp= -1-5.241e-12*x**3+2.254e-8*x**2-2.805e-5*x+1.117 + (4.e-4*x-0.1)/(math.exp(0.01*(x-trs))+1)
	return 1./(1.-B*tmp)





