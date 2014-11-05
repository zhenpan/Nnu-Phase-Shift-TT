#da#
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
import camb4py
from classy import Class

camb = camb4py.load('/home/zhenpan/camb/camb')
cosmo= Class()

def bbn(As, ns, Tau, H0, omgb, omgc, nnu):
	params = {'output': 'tCl', 'l_max_scalars':2000, 
		'N_ur':nnu, 'YHe':'BBN',
		'A_s':As, 'n_s':ns, 'k_pivot':0.05, 'tau_reio':Tau, 
		'H0':H0, 'omega_b':omgb,'omega_cdm':omgc}
	cosmo.set(params)
	cosmo.compute()
	x = cosmo.ionization_fraction(80000) 
	return (2*x-2)/(2*x-1)         #return YHe


def run(As, ns, Tau, H0, omgb, omgc, nnu, Yp):
	result = camb(**{'get_scalar_cls':True, 'temp_camb':2.7255, 
			'massless_neutrinos': nnu, 'massive_neutrinos': 0, 'helium_fraction':Yp,
			'ombh2': omgb, 'omch2': omgc, 'hubble': H0, 're_optical_depth': Tau,
			'scalar_spectral_index(1)': ns,'scalar_amp(1)': As, 'pivot_scalar':0.05})
	Out = result['misc']
	return Out


def params(ths_fid, thd_fid, nnu):
	omgb = omgb_fid
	omgc = (2.4719+0.5614*nnu)/(2.4719+0.5614*nnu_fid)*(omgb_fid+omgc_fid)-omgb
	hub0 = 50.
	hub1 = 80.
	
	while(abs(hub1/hub0 -1) > 1.e-4):
		hub = (hub0+hub1)/2.
		Yp  = bbn(As, ns, Tau, hub, omgb, omgc, nnu)
		Out = run(As, ns, Tau, hub, omgb, omgc, nnu, Yp)
		ths = float(Out['100*theta'])	
		thd = float(Out['100*theta_D'])
		zeq = float(Out[' '])	
		if(ths > ths_fid):
			hub1 = hub
		else:
			hub0 = hub
		print hub, omgb, omgc, Yp, (ths/ths_fid-1), (thd/thd_fid-1) 
	return hub, omgb, omgc, Yp, ths, thd, zeq




As  = 2.215e-9
ns  = 0.9619
Tau = 0.0925

omgb_fid= 0.02203
omgc_fid= 0.1204
nnu_fid = 3.046
H0_fid  = 67.04
Yp_fid  = bbn(As, ns, Tau, H0_fid, omgb_fid, omgc_fid, nnu_fid)

nnu1 = 1
nnu2 = 3
nnu3 = 5

Out = run(As, ns, Tau, H0_fid, omgb_fid, omgc_fid, nnu_fid, Yp_fid)
ths_fid  = float(Out['100*theta'])
thd_fid  = float(Out['100*theta_D'])	
zeq_fid  = float(Out[' '])

H01, omgb1, omgc1, Yp1, ths1, thd1, zeq1 = params(ths_fid, thd_fid, nnu1)
H02, omgb2, omgc2, Yp2, ths2, thd2, zeq2 = params(ths_fid, thd_fid, nnu2)
H03, omgb3, omgc3, Yp3, ths3, thd3, zeq3 = params(ths_fid, thd_fid, nnu3)

print '100*theta_star:', ths_fid, ths1, ths2, ths3
print '100*theta_D:', thd_fid, thd1, thd2, thd3
print 'z_eq:', zeq_fid, zeq1, zeq2, zeq3

par = np.array([[As, As, As, As], [ns, ns, ns, ns], [Tau,Tau, Tau, Tau], [H0_fid, H01, H02, H03], [omgb_fid, omgb1, omgb2, omgb3], 
[omgc_fid, omgc1, omgc2, omgc3], [nnu_fid, nnu1, nnu2, nnu3], [Yp_fid,Yp1, Yp2, Yp3]])
np.savetxt('params_bbn', par)


