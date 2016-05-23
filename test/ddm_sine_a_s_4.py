import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys

H=512
N=4096
n=np.arange(N)
Fs=16000.
f0=1000.
x=np.exp(2.*np.pi*1j*n/Fs*f0)
x+=np.random.standard_normal(N)*np.sqrt(1)
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
