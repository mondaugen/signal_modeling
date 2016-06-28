# Analyse frames using the DDM.
# Resynthesize using a modified McAulay and Quatieri method, taking into
# consideration the estimated frequency slope.
# Here we use a cubic polynomial of phase and amplitude

import numpy as np
import matplotlib.pyplot as plt
import os
import sigmod as sm

show_plots=True

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='mq_mod_cubic'

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# Hop size
H=256
# Analysis window / FFT size
N=1024
# Sample rate
Fs=16000.

# Time points freq
t_t=np.r_[0.,0.25,0.5]*1.
# signal frequency breakpoints (Hz)
f_t=np.r_[100.,200.,100.]
#f_t=np.r_[100.,102.,100.]
# Time points 
t_a=np.r_[0.,0.1,0.3,0.5]*1.
# sample indices at time points
m_t=t_t*Fs
# sample indices at time points
m_t_a=t_a*Fs
# signal amplitude breakpoints (dB)
a_t=np.r_[-10,0,0,-10]
# fit polynomial to amplitude function
d_a=np.polyfit(m_t_a,np.log(10.**(a_t/20.)),len(a_t)-1)
# Length of signal, seconds
T_x=0.5
# Length in samples
M=int(np.floor(Fs*T_x))+N
# sample indices
m=np.arange(M)
# angular velocities
w_t=f_t/Fs*2.*np.pi
#print w_t
# initial phase
phi_0=0.

# Fit polynomial to frequency function
d_=np.polyfit(m_t,w_t,len(w_t)-1)
# Phase function integral of frequency function with initial phase as smallest
# coefficient
d=np.polyint(d_,k=phi_0)
# Synthesize signal
arg_ph_x=np.polyval(d,m)
with open(plotoutpath+'_arg_ph_x.f64','w') as f:
    arg_ph_x.tofile(f)
arg_a_x=np.polyval(d_a,m)
with open(plotoutpath+'_arg_a_x.f64','w') as f:
    arg_a_x.tofile(f)
x=np.exp(1j*arg_ph_x)*np.exp(arg_a_x)

# Estimated parameters
th=[]
# Analysis window indices
n=np.arange(N)
W,dW=sm.w_dw_sum_cos(N,'c1-blackman-4')

# Plot
plt.figure(1)
plt.specgram(x,NFFT=N,noverlap=(N-H),Fs=Fs,cmap="Greys")
plt.title('Original signal: (spectrogram)')
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.gca().set_xlim(0,(len(x)-N)/float(Fs))
plt.gca().set_ylim(f_t.min()*0.5,f_t.max()*1.5)
plt.savefig(plotoutpath+'_original_spec.eps')

for h in np.arange(0,M-N,H):
    x_=x[h:h+N]
    a_=sm.ddm_p2_1_3(x_,W,dW)
    # Store estimated parameters
    th.append(a_)

# Synthesize using modified McAulay & Quatieri cubic phase method
h=0
y=np.zeros(len(x)).astype('complex_')
# Argument x of phase function exp(j*x)
arg_ph=np.zeros(len(x)).astype('double')
# Argument x of amplitude function exp(j*x)
arg_a=np.zeros(len(x)).astype('double')
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
    arg_ph[h:h+H]=np.polyval(c,np.arange(H))
    y[h:h+H]=np.exp(1j*arg_ph[h:h+H])
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
    arg_a[h:h+H]=np.polyval(d,np.arange(H))
    y[h:h+H]*=np.exp(arg_a[h:h+H])
    h+=H

# Save phase and log-amplitude polynomials
with open(plotoutpath+'_arg_ph.f64','w') as f:
    arg_ph.tofile(f)
with open(plotoutpath+'_arg_a.f64','w') as f:
    arg_a.tofile(f)

plt.figure(2)
plt.specgram(y,NFFT=N,noverlap=(N-H),Fs=Fs,cmap="Greys")
plt.title('Estimated signal (spectrogram)')
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.gca().set_xlim(0,(h-N)/float(Fs))
plt.gca().set_ylim(f_t.min()*0.5,f_t.max()*1.5)
plt.savefig(plotoutpath+'_estimated_spec.eps')

plt.figure(3)
# Plot length
N_plt_0=2000
N_plt_1=3000
plt.plot(m/float(Fs),np.real(x),c='k',label='True')
plt.plot(m/float(Fs),np.real(y),c='Gray',label='Estimated')
plt.gca().set_xlim(N_plt_0/float(Fs),N_plt_1/float(Fs))
plt.title('True vs. Estimated signal (real part)')
plt.ylabel('Amplitude')
plt.xlabel('Time (seconds)')
plt.legend()
plt.savefig(plotoutpath+'_orig_vs_est.eps')
plt.figure(4)
plt.plot(m/float(Fs),20.*np.log10(np.abs(y-x)),c='k')
plt.gca().set_xlim(0,(h-N)/float(Fs))
plt.title('Error signal (db Error)')
plt.ylabel('Amplitude (dB power)')
plt.xlabel('Time (seconds)')
plt.savefig(plotoutpath+'_error.eps')

if (show_plots):
    plt.show()
