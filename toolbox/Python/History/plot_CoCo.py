import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import datadir,outdir
plt.ion()

plt.figure()
for color,m in [('r','11'),('g','12'),('b','13')]:
  halo=np.loadtxt(datadir+'CoCo/halo'+m+'.stats')
  sub=np.loadtxt(datadir+'CoCo/sat'+m+'.stats')
  r=sub[:,1]
  f=(r<1)
  f2=(r>0.2)
  x=r[f]
  ymean=np.exp(np.interp(np.log(r), np.log(halo[:,1]), np.log(halo[:,3])))
  ymed=np.exp(np.interp(np.log(r), np.log(halo[:,1]), np.log(halo[:,8])))
  #plt.errorbar(x, sub[f,3]/ymean, sub[f, 14]/ymean, fmt='r')
  sig=sub[:, 16]
  #sig=sub[:, 18]-sub[:,16]
  #sig=sub[:, 8]/np.sqrt(sub[:,19])
  norm=(sub[f,8]/ymed[f])[-1]
  plt.errorbar(x, (sub[f,8]/ymed[f])/norm, (sig/ymed)[f]/norm, color=color, fmt='o')
  plt.fill_between(x[f], sub[f,11]/ymed[f]/norm, sub[f,12]/ymed[f]/norm, color=color, alpha=0.2)
  polypar,polyV=np.polyfit(np.log(r[f&f2]), np.log(sub[:,8]/ymed/norm)[f&f2],1, w=1/(sig/ymed)[f&f2], cov=True)
  plt.plot(x, np.exp(np.polyval(polypar, np.log(x))), '--', color=color, label=m+r':$\gamma=%.2f\pm %.2f$'%(polypar[0],np.sqrt(polyV[0,0])))
  plt.xscale('log')
  plt.yscale('log', nonposy='clip')
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$\rho_{\rm sub}/\rho_{\rm host}$')
plt.legend(loc=4)
#plt.ylim([5e-4,1e-1])
plt.plot(x, x**1.3/x[f][-1]**1.3, 'k-', label=r'AqA1: $\gamma=1.18$')


#plt.savefig(outdir+'CoCoSubprof_rat.pdf')