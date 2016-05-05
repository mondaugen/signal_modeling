import numpy as np
import matplotlib.pyplot as plt
f=400.
fs=16000.
T=1.
N_fft=1024
L=512
H=1
t=np.arange(0.,T,1./fs)
N=len(t)
n=np.mat(np.arange(0,N-H,H)).T+np.arange(H)
x=np.cos(2.*np.pi*f*t)
X=x[n]
Y=X*np.pad(np.hanning(L),(N_fft-L)/2,'constant',constant_values=0)
S=np.fft.fft(Y)
plt.figure(1)
plt.plot(t,x)
plt.figure(2)
n_plt=np.arange(0,N-H,H)
s_plt=np.angle(S[n_plt,26])
s_last=0
for s_ in np.nditer(s_plt,op_flags=['readwrite']):
    if (s_ - s_last) >= np.pi/2.:
        s_[...]-=np.pi
    if (s_ - s_last) <= -np.pi/2.:
        s_[...]+=np.pi
    s_last=s_[...]

print s_plt

plt.plot(n_plt,s_plt)
plt.show()
