from nfw import NFWHalo
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

#H=NFWHalo(183,15)
#H=NFWHalo(83,9)
H=NFWHalo(100, 7)
R=np.logspace(-2,0.3, 50)*H.Rv

msat=1e-4*H.M
cref=NFWHalo(msat).C
k=1

y=H.strip_func(NFWHalo(msat), R, k=k)
f=(R>0.01)&(R<1)
par=np.polyfit(np.log(R[f]), np.log(y)[f],1)
print par

plt.loglog(R/H.Rv, y, 'r')
#plt.plot(R/H.Rv, np.exp(np.polyval(par, np.log(R)))/h.M, '-', color=color)
plt.loglog(R/H.Rv, H.strip_func(NFWHalo(msat, cref*10**0.15), R, k=k), 'g')
plt.loglog(R/H.Rv, H.strip_func(NFWHalo(msat, cref/10**0.15), R, k=k), 'b')

plt.plot(R/H.Rv, np.exp(1.18*np.log(R/H.Rv)-0.2), 'k-', label=r'AqA1: $\gamma=1.18$')

plt.figure()
for c in [6, 15]:
  for m in [1e2, 5e4]:
	H=NFWHalo(m,c)
	R=np.logspace(-2,0.3, 50)*H.Rv
	y=H.strip_func(NFWHalo(1e-4*H.M), R, k=1)
	plt.loglog(R/H.Rv, y, label='%.0e,%d'%(m,c))
	
plt.legend(loc=4)