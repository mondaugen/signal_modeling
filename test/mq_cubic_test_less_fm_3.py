# Analyse frames using the DDM.
# Resynthesize using a modified McAulay and Quatieri method, taking into
# consideration the estimated frequency slope.
# Here we use a quartic polynomial of phase and amplitude

import numpy as np
import matplotlib.pyplot as plt
import os
import sigmod as sm
import neplot as nep

show_plots=True

# Color contrast config
# values further from 1, more contrast
clr_gamma=4.
clr_mapper=nep.PowerNormalize(clr_gamma)

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='mq_mod_quintic'

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
# Hop size
H=256
# Analysis window / FFT size
N=1024
# Analysis window indices
n=np.arange(N)
W,dW=sm.w_dw_sum_cos(N,'c1-blackman-4')

# Plot
plt.figure(1)
plt.specgram(x,NFFT=N,noverlap=(N-H),Fs=Fs,norm=clr_mapper,cmap="Greys")
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

# Polynomial error bounds
eb_c5=[]
eb_c4=[]
eb_c3=[]
eb_d5=[]
eb_d4=[]
eb_d3=[]
eb_ph=[]
eb_a=[]

# Synthesize using modified McAulay & Quatieri quintic phase method
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
    M=np.round((20.*H*(w_i0+w_i1)+(H**2.)*(psi_i0-psi_i1)+40.*(phi_i0-phi_i1))/(80.*np.pi))
    # Compute phase polynomial coefficients
    c5_p=np.r_[
            6.*(phi_i1-phi_i0+2.*np.pi*M),
            -3.*(w_i1+w_i0),
            0.5*(psi_i1-psi_i0),
            0.,
            0.,
            0]
    c4_p=np.r_[
            15.*(phi_i0-phi_i1-2.*np.pi*M),
            7.*w_i1+8.*w_i0,
            1.5*psi_i0-psi_i1,
            0.,
            0.]
    c3_p=np.r_[
            10.*(phi_i1-phi_i0+2.*np.pi*M),
            -4.*w_i1-6.*w_i0,
            0.5*psi_i1-1.5*psi_i0,
            0.]
    c5,eb_c5_=sm.polyval_mu(c5_p,1./H)
    c4,eb_c4_=sm.polyval_mu(c4_p,1./H)
    c3,eb_c3_=sm.polyval_mu(c3_p,1./H)
    eb_c5.append(eb_c5_)
    eb_c4.append(eb_c4_)
    eb_c3.append(eb_c3_)
    c2=0.5*psi_i0
    c1=w_i0
    c0=phi_i0
    c=np.r_[c5,c4,c3,c2,c1,c0]
    # evaluate phase polynomial
    ph_,eb_ph_=sm.polyval_mu(c,np.arange(H))
    arg_ph[h:h+H]=ph_
    eb_ph += list(eb_ph_)
    y[h:h+H]=np.exp(1j*ph_)
    # compute amplitude polynomial coefficients
    a0_i0=np.real(th[i][0])
    a0_i1=np.real(th[i+1][0])
    a1_i0=np.real(th[i][1])
    a1_i1=np.real(th[i+1][1])
    a2_i0=np.real(th[i][2])
    a2_i1=np.real(th[i+1][2])
    # Find coefficients of quintic amplitude function:
    # mu(t) = d5*t^5 + d4*t^4 + d3*t^3 + d2*t^2 + d1*t + d0
    # at t=0, mu(0)= a0_i0 so
    d0=a0_i0
    d5_p=np.r_[
            6.*(a0_i1-a0_i0),
            -3.*(a1_i1+a1_i0),
            0.5*(a2_i1-a2_i0),
            0.,
            0.,
            0.]
    d4_p=np.r_[
            15.*(a0_i0-a0_i1),
            7.*a1_i1+8.*a1_i0,
            1.5*a2_i0-a2_i1,
            0.,
            0.]
    d3_p=np.r_[
            10.*(a0_i1-a0_i0),
            -4.*a1_i1-6.*a1_i0,
            0.5*a2_i1-1.5*a2_i0,
            0.]
    d5,eb_d5_=sm.polyval_mu(d5_p,1./H)
    d4,eb_d4_=sm.polyval_mu(d4_p,1./H)
    d3,eb_d3_=sm.polyval_mu(d3_p,1./H)
    eb_d5.append(eb_d5_)
    eb_d4.append(eb_d4_)
    eb_d3.append(eb_d3_)
    d5=np.polyval(d5_p,1./H)
    d4=np.polyval(d4_p,1./H)
    d3=np.polyval(d3_p,1./H)
    d2=0.5*a2_i0
    d1=a1_i0
    d0=a0_i0
    d=np.r_[d5,d4,d3,d2,d1,d0]
   # Multiply by amplitude function
    arg_a[h:h+H],eb_a_=sm.polyval_mu(d,np.arange(H))
    y[h:h+H]*=np.exp(arg_a[h:h+H])
    eb_a += list(eb_a_)
    h+=H

# Save phase and log-amplitude polynomials
with open(plotoutpath+'_arg_ph.f64','w') as f:
    arg_ph.tofile(f)
with open(plotoutpath+'_arg_a.f64','w') as f:
    arg_a.tofile(f)

# Save true and estimated signals
with open(plotoutpath+'_true_x.dat','w') as f:
    x.tofile(f)
with open(plotoutpath+'_est_x.dat','w') as f:
    y.tofile(f)


plt.figure(2)
plt.specgram(y,NFFT=N,noverlap=(N-H),Fs=Fs,norm=clr_mapper,cmap="Greys")
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
plt.title('Original vs. estimated signal: error')
plt.ylabel('Amplitude (dB power)')
plt.xlabel('Time (seconds)')
plt.savefig(plotoutpath+'_error.eps')
plt.figure(5)
tmp=np.array([np.log(eb_ph_)/np.log(10.) for eb_ph_ in
    eb_ph])
ma_,mai_=sm.lextrem(tmp,comp='max')
plt.plot(np.arange(len(eb_ph))[mai_],tmp[mai_],label="Phase",c='k',ls='-')
tmp=np.array([np.log(eb_a_)/np.log(10.) for eb_a_ in eb_a])
ma_,mai_=sm.lextrem(tmp,comp='max')
plt.plot(np.arange(len(eb_a))[mai_],tmp[mai_],
        label="Amplitude",c='k',ls=':')
plt.ylabel('Absolute error bound ($\log_{10}$)')
plt.title('Polynomial evaluation error bound')
plt.legend(loc='best')
plt.savefig(plotoutpath+'_poly_eval_err.eps')
plt.figure(6)
plt.plot(np.arange(len(eb_c5)),np.log(np.array([eb_c5,eb_c4,eb_c3,eb_d5,eb_d4,eb_d3]).T)/np.log(10.))
plt.xlabel('Frame number')
plt.ylabel('Absolute error bound')
plt.title('Polynomial evaluation error bound')

if (show_plots):
    plt.show()
