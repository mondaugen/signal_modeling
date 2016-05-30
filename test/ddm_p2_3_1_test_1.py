import numpy as np
import sigmod as sm
a0=0
a1=complex('1e-3+0.5j')
a2=complex('-1e-6+1e-4j')
N=512
H=256
L=2048
l=np.arange(L)
x=np.exp(0+a1*l+a2*l**2.)
n=np.arange(N)
w=0.5-0.5*np.cos(2.*np.pi*n/N)
dw=np.pi/N*np.sin(2.*np.pi*n/N)
a=sm.ddm_p2_3_1(x,H,w,dw,2)
print a

