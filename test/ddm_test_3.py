import numpy as np
import matplotlib.pyplot as plt

j1=complex('1j')
a0=0
a1=complex('1e-3+0.50j')
N=513
n=np.arange(-(N-1)/2,(N-1)/2+1)
x=np.exp(a0+a1*n)
x+=np.random.standard_normal(N)*np.sqrt(0)
h=0.5+0.5*np.cos(2.*np.pi*n/N)
dh=-np.pi/N*np.sin(2.*np.pi*n/N)
F=np.exp(np.outer(n,j1*np.pi*-2.*n/N))
Fh=F*h
Fdh=F*(np.outer(-2.*np.pi*j1*n/N,h)+dh)
#Fdh=F*dh
Xh=np.inner(Fh,x)
Xdh=np.inner(Fdh,x)
a=-Xdh/Xh
print a[np.abs(Xh).argmax()]
plt.figure(1)
plt.plot(n,np.real(x))
plt.figure(2)
plt.plot(n,np.abs(Xh))
plt.show()
