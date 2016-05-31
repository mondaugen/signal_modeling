# Compare 3 frames 1 bin, 1 frame 3 bins and 3 frames 3 bins
import numpy as np
import sigmod as sm
a0=0
a1=complex('0+0j')
a2=complex('-1e-6+1e-4j')
#a2=complex('0+1e-4j')
N=512
H=128
L=2048
l=np.arange(L)-H
x=np.exp(a0+a1*l+a2*l**2.)
n=np.arange(N)
w=0.5-0.5*np.cos(2.*np.pi*n/N)
dw=np.pi/N*np.sin(2.*np.pi*n/N)
a=sm.ddm_p2_3_1(x,H,w,dw,2)
print a
a_=sm.ddm_p2_1_3(x[H:],w,dw)
print a_
a__=sm.ddm_p2_3_3(x,H,w,dw,2)
print a__

