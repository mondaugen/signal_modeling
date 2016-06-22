# Take DFT of signal in frames and then resynthesize by interpolating using the
# McAulay and Quatieri method.
import numpy as np
import matplotlib.pyplot as plt
import os

show_plots=True

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='mq_cubic'

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# Hop size
H=256
# Analysis window / FFT size
N=1024
# Time points
t_t=np.r_[0.,0.25,0.5]*1.
# signal starts at 100 Hz, goes to 110Hz and ends at 105Hz
f_t=np.r_[100.,500.,200.]
# Amplitude constant
A_t=np.ones(len(f_t))
# Sample rate
Fs=16000.
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
# sample indices at time points
m_t=t_t*Fs

# Fit polynomial to frequency function
d_=np.polyfit(m_t,w_t,len(w_t)-1)
# Phase function integral of frequency function with initial phase as smallest
# coefficient
d=np.polyint(d_,k=phi_0)
# Synthesize signal
x=np.exp(1j*np.polyval(d,m))

# Estimated parameters
th=[]
# Analysis window
W=np.hanning(N+1)[:N]
W_0=np.sum(W)
# Total length of FFT (for zero padding)
N_fft=2*N

# Plot
plt.figure(1)
plt.specgram(x,NFFT=N,noverlap=(N-H),Fs=Fs,cmap="Greys")
plt.title('Original signal: (spectrogram)')
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.gca().set_xlim(0,(len(x)-2*N)/float(Fs))
plt.gca().set_ylim(f_t.min()*0.5,f_t.max()*2.)
plt.savefig(plotoutpath+'_original_spec.eps')

for h in np.arange(0,M-N,H):
    x_=x[h:h+N]*W
    X_=np.fft.fft(x_,N_fft)/W_0
    kma=np.abs(X_).argmax()
    if (kma > 0) and (kma < (N_fft-1)):
        # Estimate better maximum using quadratic interpolation
        B_X=np.matrix([
            [0.5,-1.,0.5],
            [-0.5,0,0.5],
            [0,1,0]
        ])
        b_X=np.matrix([
            [np.abs(X_)[kma-1]],
            [np.abs(X_)[kma]],
            [np.abs(X_)[kma+1]]
        ])
        a_X=B_X*b_X
        a_X=np.array(a_X).flatten()
        kma_=-a_X[1]/(2.*a_X[0])
        Xma_=np.polyval(a_X,kma_)
        kma+=kma_
    else:
        Xma_=np.abs(X_)[kma]
    # Store phase, frequency, amplitude
    th.append([np.angle(X_[kma]),kma/float(N_fft)*2.*np.pi,Xma_])

# Synthesize using McAulay & Quatieri cubic phase method
h=0
y=np.zeros(len(x)).astype('complex_')
for i in xrange(len(th)-1):
    phi_i0=th[i][0]
    phi_i1=th[i+1][0]
    w_i0=th[i][1]
    w_i1=th[i+1][1]
    A_i0=th[i][2]
    A_i1=th[i+1][2]
    # Compute M*
    M=np.round((H/2.*(w_i1-w_i0)-(phi_i1-phi_i0-w_i0*H))/(2.*np.pi))
    # Compute polynomial coefficients
    q_=phi_i1-phi_i0-w_i0*H+2.*np.pi*M
    r_=w_i1-w_i0
    c=np.r_[
            -2.*q_/(H**3) + r_/(H**2),
            3.*q_/(H**2) - r_/H,
            w_i0,
            phi_i0
    ]
    y[h:h+H]=np.exp(1j*np.polyval(c,np.arange(H)))
    # Linearly interpolate amplitude
    y[h:h+H]*=np.interp(np.arange(H),np.r_[0,H],np.r_[A_i0,A_i1])
    h+=H

plt.figure(2)
plt.specgram(y,NFFT=N,noverlap=(N-H),Fs=Fs,cmap="Greys")
plt.title('Estimated signal (spectrogram)')
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.gca().set_xlim(0,(h-N)/float(Fs))
plt.gca().set_ylim(f_t.min()*0.5,f_t.max()*2.)
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
