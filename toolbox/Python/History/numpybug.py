import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
sample_size=1e5
xmin,xmax=-10.,0.
data=np.random.rand(sample_size)*(xmax-xmin)+xmin #some uniform data
weight=np.exp(-4*data) #weights

y,x=np.histogram(data,50,weights=weight)
print y[-1]
print np.sum(weight[(data>x[-2])&(data<=x[-1])])
#yy=st.histogram(data,50,weights=weight)[0]
#print yy[-1]

bin_index=np.digitize(data, x)
yy=np.bincount(bin_index, weights=weight)[1:-1]
print yy[-1]