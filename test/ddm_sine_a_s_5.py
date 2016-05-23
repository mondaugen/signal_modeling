import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys

H=512
N=4096
n=np.arange(N)
Fs=16000.
# start frequency
f0=1000.
# goal frequency
f1=1100.
# goal time
N1=1024.
psi=(f1-f0)/N1*2.*np.pi/Fs
x=np.exp(1j*(2.*np.pi*n/Fs*f0+0.5*psi*n**2.))
x+=np.random.standard_normal(N)*np.sqrt(.1)
y=np.zeros(N).astype('d')
#plt.plot(n,np.real(x))
a=[]
h=0
M=2048
m=np.arange(M)
w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
while ((h+M) <= N):
    a.append(sm.ddm_p2_1_3(x[h:(h+M)],w,dw))
    print np.imag(a[-1][1])*Fs/(2.*np.pi)
    h+=H
