# Take DFT of signal in frames and then resynthesize by interpolating using the
# McAulay and Quatieri method.
import numpy as np
import matplotlib.pyplot as plt

# Time points
t_t=np.r_[0.,0.25,0.5]*10.
# signal starts at 100 Hz, goes to 110Hz and ends at 105Hz
f_t=np.r_[100.,300.,200.]
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
print w_t
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
H=128
# Analysis window / FFT size
N=512
# Analysis window
W=np.hanning(N+1)[:N]
# Total length of FFT (for zero padding)
N_fft=N
for h in np.arange(0,M-N,H):
    x_=x[h:h+N]*W
    X_=np.fft.fft(x_,N_fft)
    kma=np.abs(X_).argmax()
    # Store phase, frequency, amplitude
    th.append([np.angle(X_[kma]),kma/float(N_fft)*2.*np.pi,np.abs(X_[kma])])

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
plt.specgram(y,Fs=Fs)
plt.title('Estimated signal: spectrogram')

plt.show()
