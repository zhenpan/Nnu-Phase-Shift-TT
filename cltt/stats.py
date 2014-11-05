
import math
import numpy as np
import matplotlib.pyplot as plt

err  = np.loadtxt('res1.dat')
mean = np.zeros(18)
var  = np.zeros(18)
count= np.zeros(18)

for i in range(18):
	res  = [x for x in err[:][i]]
	Num  = len(res)
	mean[i] = sum(res)/Num
	res2 = [x**2-mean[i]**2 for x in res]
	var2 = sum(res2)/Num
	var[i]  = math.sqrt(var2)
	count[i] = i+1

#print mean, var
plt.figure(0)
plt.errorbar(count, mean, yerr = var, fmt='o')
plt.plot((0,18), (0,0))
plt.xlabel(r'$i\rm th \ peak\ or\ valley$')
plt.ylabel(r'${\rm residual}\  \delta\ell$')
plt.show()


mean = np.zeros(100)
var  = np.zeros(100)
count= np.zeros(100)

for i in range(100):
	res = [err[i][j] for j in range(12)]
	Num = len(res)
	mean[i] = sum(res)/Num
	res2 = [x**2-mean[i]**2 for x in res]
	var2 = sum(res2)/Num
	var[i]= math.sqrt(var2)
	count[i] = i+1	

mmean = sum(mean)/100
tmp = [x**2-mmean**2 for x in mean]
vvar = math.sqrt(sum(tmp)/100)
print mmean, vvar

plt.figure(1)
plt.plot(count, mean, 'o')
#plt.errorbar(count, mean, yerr = var, fmt='o')
plt.plot((0,100),(0,0))
plt.xlabel(r'$i\rm th \ pairs \ of \ model$')
plt.ylabel(r'${\delta\bar\ell: \rm average\ of \ residuals \ for \ \ell < 2000}$')
plt.title(r'$\delta\bar\ell = 0.01\pm 0.2$')
plt.show()
	
