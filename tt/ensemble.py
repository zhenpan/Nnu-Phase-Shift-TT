#da#
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
import random
import camb4py
from cosmoslik import chains

camb = camb4py.load('/home/zhenpan/camb/camb')

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
		Yp  = Helium(thd_fid, omgb, omgc, nnu, hub)
		Out = run(As, ns, Tau, hub, omgb, omgc, nnu, Yp)
		ths = float(Out['100*theta'])	
		thd = float(Out['100*theta_D'])
		zeq = float(Out[' '])	
		if(ths > ths_fid):
			hub1 = hub
		else:
			hub0 = hub
		#print hub, omgb, omgc, Yp, (ths/ths_fid-1), (thd/thd_fid-1) 
	return hub, omgb, omgc, Yp, ths, thd, zeq

def Helium(thd_fid, omgb, omgc, nnu, H0):
	Y0 = 0.1
	Y1 = 0.4
	thd= 2*thd_fid

	while(abs(thd/thd_fid-1) > 1.e-4 and abs(Y1/Y0-1) > 1.e-5 ):
		Yp  = (Y0+Y1)/2.
	 	Out = run(As, ns, Tau, H0, omgb, omgc, nnu, Yp)
		thd = float(Out['100*theta_D'])	
		if(thd > thd_fid):
			Y1 = Yp
		else:
			Y0 = Yp
	return Yp 


data = chains.load_chain('data_thin100')

Amp    = data['A*']
Ndex   = data['ns']
T_reio = data['tau']
omgbh2 = data['omegabh2']
omgch2 = data['omegach2']
Hubble = data['H0*']

Nsamp= 100
par1 = np.zeros((Nsamp,8))
par2 = np.zeros((Nsamp,8))

for i in range(Nsamp):
	k   = 10*i
	nnu1 = random.randrange(1,4) #[1 2 3]
	nnu2 = nnu1+3		     #[]+3	
	#nnu1 = random.randrange(1,6) #random[1 2 3 4 5]
	#nnu2 = random.randrange(1,6)+0.5 #avoiding same value

	As  = Amp[k]*1e-9
	ns  = Ndex[k]
	Tau = T_reio[k] 
	omgb_fid= omgbh2[k]
	omgc_fid= omgch2[k]
	H0_fid  = Hubble[k] 
	nnu_fid = 3.046
	Yp_fid  = 0.2477

	Out = run(As, ns, Tau, H0_fid, omgb_fid, omgc_fid, nnu_fid, Yp_fid)
	ths_fid  = float(Out['100*theta'])
	thd_fid  = float(Out['100*theta_D'])	
	zeq_fid  = float(Out[' '])

	#ths_fid2  = data['thetastar*'][k]
	#thd_fid2  = data['thetad*'][k]
	#zeq_fid2  = data['zeq*'][k]
	#print '100*theta_star:', ths_fid, ths_fid2
	#print '100*theta_D:', thd_fid, thd_fid2
	#print 'z_eq:', zeq_fid, zeq_fid2

	H01, omgb1, omgc1, Yp1, ths1, thd1, zeq1 = params(ths_fid, thd_fid, nnu1)
	H02, omgb2, omgc2, Yp2, ths2, thd2, zeq2 = params(ths_fid, thd_fid, nnu2)
	
	print 'i_th:', i
	print 'Nnu:', nnu_fid, nnu1, nnu2
	print '100*theta_star:', ths_fid, ths1, ths2
	print '100*theta_D:', thd_fid, thd1, thd2
	print 'z_eq:', zeq_fid, zeq1, zeq2
	
	params1 = [As, ns, Tau, H01, omgb1, omgc1, nnu1, Yp1]
	params2 = [As, ns, Tau, H02, omgb2, omgc2, nnu2, Yp2]
	for j in range(8):
		par1[i][j] = params1[j]			
		par2[i][j] = params2[j]			

np.savetxt('params3a', par1)
np.savetxt('params3b', par2)


