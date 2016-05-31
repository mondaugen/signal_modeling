# The following shows that if a decent estimate for the centre frequencies of
# two close sinusoids can be made, then the function of their phases has a local
# minimum that corresponds to a good set of initial phases.

import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

N=4000
n=np.arange(N)
w1=0.01
w2=0.0121
phi1=2
phi2=1.2
x1=np.cos(phi1+2*np.pi*w1*n)
x2=0.5*np.cos(phi2+2*np.pi*w2*n)
x=x1+x2+np.random.standard_normal(N)*np.sqrt(1)
plt.figure(1)
plt.plot(n,x1,n,x2,n,x)

phi1_=np.linspace(-np.pi,np.pi,100)
phi2_=np.linspace(-np.pi,np.pi,100)

R=np.ndarray((len(phi1_),len(phi2_)))
for i in xrange(len(phi1_)):
    for j in xrange(len(phi2_)):
        R[i,j]=np.sum((x-np.cos(phi1_[i]+2.*np.pi*w1*n)
            -np.cos(phi2_[j]+2.*np.pi*w2*n))**2.)

plt.figure(2)
plt.imshow(R,extent=[phi2_.min(),phi2_.max(),phi1_.min(),phi1_.max()])
plt.title('Mean square error minimizing grid-search')

w1_=np.linspace(0,np.pi,100)
w2_=np.linspace(0,np.pi,100)
R2=np.ndarray((len(w1_),len(w2_)))
for i in xrange(len(w1_)):
    for j in xrange(len(w2_)):
        R2[i,j]=np.sum((x-np.cos(phi1+w1_[i]*n)
            -np.cos(phi2+w2_[j]*n))**2.)

plt.figure(3)
plt.imshow(R2)

plt.show()
