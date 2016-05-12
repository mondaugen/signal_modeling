import numpy as np
import matplotlib.pyplot as plt
N=1024.
n=np.arange(N+1)
Fs=16000.
f=302.
x=np.cos(n/Fs*2.*np.pi*f)
x+=np.random.standard_normal(N+1)*np.sqrt(1e0)
X_0=np.fft.fft(x[:N]*np.hanning(N))
X_1=np.fft.fft(x[1:]*np.hanning(N))
n_bin=np.round(f/(Fs/N))
print np.angle(X_1[n_bin]/X_0[n_bin])/(np.pi*2.)*Fs
