# Analyse frames using the DDM.
# Resynthesize using a modified McAulay and Quatieri method, taking into
# consideration the estimated frequency slope.

import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

# Time points
t_t=np.r_[0.,0.25,0.5]*10.
# signal starts at 100 Hz, goes to 110Hz and ends at 105Hz
f_t=np.r_[100.,1000.,200.]
# Amplitude constant
A_t=np.ones(len(f_t))
# Sample rate
Fs=16000.
# Length of signal, seconds
T_x=10.
# Length in samples
M=int(np.floor(Fs*T_x))
# sample indices
m=np.arange(M)
# angular velocities
w_t=f_t/Fs*2.*np.pi
#print w_t
# initial phase
phi_0=0.
# sample indices at time points
m_t=t_t*Fs

# Fit polynomial to frequency function
d_=np.polyfit(m_t,w_t,len(w_t)-1)
# Phase function integral of frequency function with initial phase as smallest
# coefficient
d=np.polyint(d_,k=phi_0)
# Synthesize signal
x=np.exp(1j*np.polyval(d,m))

# Plot
plt.figure(1)
plt.specgram(x,Fs=Fs)
plt.title('Original signal: spectrogram')

# Estimated parameters
th=[]
# Hop size
H=512
# Analysis window / FFT size
N=2048
# Analysis window indices
n=np.arange(N)
# Analysis window
wc=np.r_[0.358735,0.488305,0.141265,0.011695]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./N*np.outer(np.arange(len(wc)),n)))
W=np.sum(w_,0)
dw_=((2.*np.pi/N*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
        *np.sin(np.pi*2./N*np.outer(np.arange(1,len(wc)),n)))
dW=np.sum(dw_,0)
for h in np.arange(0,M-N,H):
    x_=x[h:h+N]
    a_=sm.ddm_p2_1_3(x_,W,dW)
    # Store estimated parameters
    th.append(a_)

# Synthesize using modified McAulay & Quatieri cubic phase method
h=0
y=np.zeros(len(x)).astype('complex_')
for i in xrange(len(th)-1):
    phi_i0=np.imag(th[i][0])
    phi_i1=np.imag(th[i+1][0])
    w_i0=np.imag(th[i][1])
    w_i1=np.imag(th[i+1][1])
    psi_i0=np.imag(th[i][2])
    psi_i1=np.imag(th[i+1][2])
    # Compute M*
    M=np.round(((w_i1+w_i0)*H-2.*(phi_i1-phi_i0))/(4.*np.pi))
    # Compute polynomial coefficients
    q_=phi_i1-phi_i0+2.*np.pi*M
    r_=0.5*(w_i1+w_i0)
    s_=0.5*(psi_i1+psi_i0)
    # Done this way for more numerical stability, maybe?
    b=np.asmatrix(np.c_[q_,r_,s_].T)
    V=np.asmatrix(np.array([
        [4./H,-4.,0],
        [-6./H,6.,0.5*H],
        [3./H,-2.,-0.5*H]
    ]))
    D=np.diag([1./(H**2.),1./H,1.])
    c_=D*V*b
    c=np.array(np.vstack((c_,phi_i0))).flatten()
    y[h:h+H]=np.exp(1j*np.polyval(c,np.arange(H)))
    a0_i0=np.real(th[i][0])
    a0_i1=np.real(th[i+1][0])
    a1_i0=np.real(th[i][1])
    a1_i1=np.real(th[i+1][1])
    a2_i0=np.real(th[i][2])
    a2_i1=np.real(th[i+1][2])
    # Find coefficients of cubic amplitude function:
    # mu(t) = d3*t^3 + d2*t^2 + d1*t + d0
    # at t=0, mu(0)= a0_i0 so
    d0=a0_i0
    # We find the rest of the coefficients using least squares
    # constrained so that
    # mu(H) = a0_i1
    # The least squares minimizes the error fit of the following:
    # mu'(0) = a1_i0
    # mu'(H) = a1_i1
    # mu''(0) = a2_i0
    # mu''(H) = a2_i1
    # Therefore the observation matrix is
    A_o=np.array([
        [0.,0.,1.],
        [3.*H*H,2*H,1.],
        [0.,2.,0.],
        [6.*H,2.,0.]
    ])
    # The constraint matrix is
    A_c=np.array([[H**3.,H**2.,H]])
    # Constraint vector
    b_c=np.array([[a0_i1-a0_i0]])
    # observation vector
    b_o=np.c_[a1_i0,a1_i1,a2_i0,a1_i1].T
    d_=sm.lstsq_c(A_o,b_o,A_c,b_c)
    d=np.array(np.vstack((d_,d0))).flatten()
    # Multiply by amplitude function
    y[h:h+H]*=np.exp(np.polyval(d,np.arange(H)))
#    y[h:h+H]*=np.interp(np.arange(H),np.r_[0,H],np.r_[a0_i0,a0_i1])
#    print 'phase error: %f' % (np.polyval(c,H)-2.*np.pi*M-phi_i1,)
#    print 'amplitude error: %f' % (np.polyval(d,H)-a0_i1,)
    h+=H

plt.figure(2)
plt.specgram(y,Fs=Fs)
plt.title('Estimated signal: spectrogram')

plt.figure(3)
plt.plot(m,np.real(x),c='b')
plt.title('True signal (real part)')
plt.plot(m,np.real(y),c='g')
plt.figure(4)
plt.plot(m,np.real(y))
plt.title('Estimated signal (real part)')
plt.figure(5)
plt.plot(m,20.*np.log10(np.abs(y-x)))
plt.title('Error signal (db Error)')

plt.show()
