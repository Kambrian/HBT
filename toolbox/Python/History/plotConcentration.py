import numpy as np
import matplotlib.pyplot as plt
from nfw import mean_concentration

plt.ion()

m=np.logspace(-6,16,20)
plt.figure()
plt.loglog(m, mean_concentration(m/1e10, 'Ludlow'), 'r-', label='Ludlow14')
plt.loglog(m, mean_concentration(m/1e10, 'MaccioW1'), 'g--', label='Maccio08')
c=mean_concentration(m/1e10, 'MaccioW1')
plt.fill_between(m, c*np.exp(0.3), c/np.exp(0.3), color='g', alpha=0.2)
plt.xlabel(r'$M_{200}[M_{\odot}/h]$')
plt.ylabel(r'$c$')
plt.legend(loc=1)
plt.savefig('/work/Projects/SubProf/plots/Concentration.pdf')	