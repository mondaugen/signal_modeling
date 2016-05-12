import numpy as np
import matplotlib.pyplot as plt
N=1024.
n=np.arange(-1,N+1)
Fs=16000.
T=1./Fs
phi=0
w=2.*np.pi*0.1
p_f_delta=0.01 # percentage of change over 1 second
psi=(w*(1.+p_f_delta)-w)
t=n/Fs
x=np.exp(1j*(phi+w*t+0.5*psi*(t**2.)))
x+=np.random.standard_normal(N+2)*np.sqrt(0)
X__T=np.fft.fft(x[:N]*np.hanning(N))
X_0=np.fft.fft(x[1:(N+1)]*np.hanning(N))
X_T=np.fft.fft(x[2:(N+2)]*np.hanning(N))
d_T_0=X_0/X__T
d_0_T=X_T/X_0
phi_0=np.angle(X_0)
phi__T=phi_0-np.angle(d_T_0)
phi_T=phi_0+np.angle(d_0_T)
A=np.array([[1.,-T,0.5*(T**2.)],
    [1.,0.,0.],
    [1.,T,0.5*(T**2.)]])
A_1=np.mat(np.linalg.inv(A))
b=A_1*np.mat(np.vstack((phi__T,phi_0,phi_T)))
print 'Original (phi,w,psi)'
print '(%f,%f,%f)' % (phi,w,psi)
print 'bin'
print w/(2.*np.pi/N)
n_bin=np.round(w/(2.*np.pi/N))
print 'Estimated'
print b.T[n_bin-1,:]
print b.T[n_bin,:]
print b.T[n_bin+1,:]
print 'Error'
print (b.T-np.mat([[phi,w,psi]]))[n_bin-1,:]
print (b.T-np.mat([[phi,w,psi]]))[n_bin,:]
print (b.T-np.mat([[phi,w,psi]]))[n_bin+1,:]
