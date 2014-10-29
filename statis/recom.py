#a#
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
from classy import Class

cosmo= Class()


def run(As, ns, Tau, H0, omgb, omgc, nnu, Yp, k):
	params = {'output': 'tCl, mPk', 'P_k_max_1/Mpc':0.2, 
		'k_output_values': k, 'YHe': Yp,
		'A_s': As, 'n_s': ns, 'k_pivot':0.05, 'tau_reio': Tau,
		'H0': H0, 'omega_b': omgb, 'omega_cdm': omgc, 'N_ur': nnu}

	cosmo.set(params)
	cosmo.compute()
        cosmo.struct_cleanup()
	return 0

def rec(a):
	i = 0
	while(a[i] < a_rec):
		i = i+1
	return i

def rs(x,omgb, omgc,nnu, H0):
	omgg  = 2.474e-5    	       
	omgr  = (2.474 + 0.562*nnu)*1.e-5    	       
	omgm  = omgb + omgc
	omgl  = (H0/100.)**2 - omgm - omgr
	ratio = 0.75*(omgb/omgg)*x	
	return 3000./math.sqrt(omgr + omgm*x + omgl*x**4)/math.sqrt(3*(1+ratio))

def recomb(As, ns, tau, H0, omgb, omgc, nnu, Yp, k, r_s):
	Cls  = run(As, ns, tau, H0, omgb, omgc, nnu, Yp, k) 
	eta, a, delta_g, psi = np.loadtxt('output/perturbations_k0_s.dat', usecols=(0,1,2,11)).T

	file = open('output/perturbations_k0_s.dat')
	line = file.readline()
	k_s  = float(line[35:54])                                           #Mpc^-1
	phi  = r_s*k_s/math.pi                                           #phase

	i_rec = rec(a)                                                    #position of recombination 
	xnew  = [a[i_rec-1+i] for i in range(4)] 
	ynew  = [0.25*delta_g[i_rec-1+i] for i in range(4)]
	#ynew  = [0.25*delta_g[i_rec-1+i] - psi[i_rec-1+i] for i in range(4)]
	f     = interp1d(xnew, ynew, kind = 'cubic')

	if(i_rec == 0):
		Delta = Delta2
	else:
		Delta = f(a_rec)	
	return k_s, phi, Delta


def samp_k(r_s, Npeak, Nsamp):  #kr_s/pi = (0.5, 1) +i
	const  = math.pi/r_s
	size = Npeak*Nsamp
	k_samp = np.zeros(size)
	for i in range(Npeak):
		for j in range(Nsamp):
			n = i*Nsamp+j
			k_samp[n] = (i + 0.5+0.5/Nsamp*j)*const 
	return k_samp
				


def shift(m):
	r_s1  = quad(rs, 0., a_rec, args=(omgbh21[m], omgch21[m], Nnu1[m], Hubble1[m]))[0]        #sound horizon at recombination 
	r_s2  = quad(rs, 0., a_rec, args=(omgbh22[m], omgch22[m], Nnu2[m], Hubble2[m]))[0]        #sound horizon at recombination 


	k_samp1 = samp_k(r_s1, Npeak, Nsamp)  
	k_samp2 = samp_k(r_s2, Npeak, Nsamp)  
	Num     = k_samp1.size

	for i in range(Num):
		j = m*Num + i
		kmode1[j], phi1[j], Delta1[j] = recomb(Amp1[m], Ns1[m], T_reio1[m], Hubble1[m], omgbh21[m], omgch21[m], Nnu1[m], Helium1[m], k_samp1[i], r_s1)
		kmode2[j], phi2[j], Delta2[j] = recomb(Amp2[m], Ns2[m], T_reio2[m], Hubble2[m], omgbh22[m], omgch22[m], Nnu2[m], Helium2[m], k_samp2[i], r_s2)
		print m, i, j, phi1[j], phi2[j]

	out1 = np.column_stack((kmode1, phi1, Delta1))
	out2 = np.column_stack((kmode2, phi2, Delta2))
	np.savetxt('data1', out1)
	np.savetxt('data2', out2)	




a_rec = 1./1000.
Amp1, Ns1, T_reio1, Hubble1, omgbh21, omgch21, Nnu1, Helium1 = np.loadtxt('params1a').T
Amp2, Ns2, T_reio2, Hubble2, omgbh22, omgch22, Nnu2, Helium2 = np.loadtxt('params1b').T

Npeak = 8   # num of peaks covered
Nsamp = 25  # num of points sampled per peak
Ncosm = 100 # num of cosmolgy pairs

Ntot = Npeak*Nsamp*Ncosm

kmode1 = np.zeros(Ntot)
kmode2 = np.zeros(Ntot)
phi1   = np.zeros(Ntot)
phi2   = np.zeros(Ntot)
Delta1 = np.zeros(Ntot)
Delta2 = np.zeros(Ntot)

for m in range(len(Amp1)):
	shift(m)
