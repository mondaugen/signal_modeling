# A test computing the 2nd order argument parameters using DDM.
import numpy as np
import matplotlib.pyplot as plt

def lmaxw(x,b,o,th):
    """
    Find local maxima by looking in windowed sections of length b and
    overlap o of array x. Only return indices whose underlying value is greater
    than th.
    """
    N=x.size
    ci=np.arange(0,N-b,o)
    ri=np.arange(b)
    S=np.add.outer(ci,ri)
    mi=np.argmax(x[S],1)
    mi=np.unique(mi+ci)
    mi=mi[np.nonzero(np.greater(x[mi],th))]
    return mi

# Synthesize test signal
# For this file, the coefficients of the modulating polynomial are with respect
# to the sample number.

# Number of partials
P=5

# time in seconds
T=10.

# Sample rate
Fs=16000

# Length in samples
N=Fs*T

# Sample indices
n=np.arange(N)

# Polynomial coefficients

# real parts
a_r=np.array([
    # initial amplitude, small range of values
    np.random.uniform(1.e-1,1.,P),
    # 1st amplitude coefficient
    # picked so that amplitude does not go up more than by a factor of 2
    # and down by a factor of 100 from the beginning to the end of the signal
    np.random.uniform(np.log(1.e-2),np.log(2.),P)/N,
    # The same but time is squared so divided by time squared
    np.random.uniform(np.log(1.e-2),np.log(2.),P)/(N**2.)
    ])
a_r=np.flipud(a_r)

# imaginary parts
a_i=np.array([
    # initial phase, between -pi and pi
    np.random.uniform(-1.*np.pi,1.*np.pi,P),
    # 1st frequency coefficient. Frequency between 0 and pi
    np.random.uniform(0.,np.pi,P),
#    np.pi*2.*10./Fs*np.ones(P),
    # 2nd frequency coefficient. Frequency doesn't change by more than octave
    # throughout the sound
    0.5*np.random.uniform(-1.,1.,P)/(float(N))
    ])
a_i[2,:]*=a_i[1,:]
a_i=np.flipud(a_i)

x=np.exp(np.apply_along_axis(np.polyval,0,a_r,n)
        + complex('1j')*np.apply_along_axis(np.polyval,0,a_i,n))
xs=np.sum(x,1)

# Add noise
xs+=np.random.standard_normal(N)

# Size of DFT in STFT
N_stft=512

# Hop size for STFT
H_stft=128

# Window in STFT
W_stft=np.hanning(N_stft)

# frame up the signal to perform STFT
n_stft=np.add.outer(np.arange(0,N-N_stft,H_stft),np.arange(N_stft))
X_stft=xs[n_stft.astype('i')]
X_stft*=W_stft
X_stft=np.fft.fft(X_stft)/np.sum(W_stft)

lmaxs=[]
maxis=[]
for row_ in np.abs(X_stft):
    z_=np.nonzero((row_[1:] > row_[:-1]) & (row_[:-1] > row_[1:]))
    maxis.append(z_)
    lmaxs.append(lmaxw(row_[z_],32,16,1e-3))

plt.figure(1)
plt.plot(n,np.abs(x))

plt.figure(2)
plt.imshow(np.abs(X_stft).T,origin='lower',interpolation='none')
for l in xrange(len(lmaxs)):
    plt.scatter([l for _ in xrange(len(lmaxs[l]))],maxis[l][lmaxs[l]],c='r')

plt.show()

