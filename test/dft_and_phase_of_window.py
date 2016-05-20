# Examine the contribution of phase by the fourier window
import numpy as np

N=32
k=2
n=np.arange(N)
x=np.exp(1j*1.+1j*2.*np.pi*n/N*k)
w1=0.5+0.5*np.cos(2*np.pi*np.arange(-N/2,N/2)/N)
X1=np.fft.fft(x*w1)
w2=0.5+0.5*np.cos(2*np.pi*np.arange(-N/2+1,N/2)/(N-1))
w2=np.r_[w2,0]
X2=np.fft.fft(x*w2)
w3=0.5+0.5*np.cos(2*np.pi*np.arange(-(N-2)/2,(N-2)/2)/(N-2))
w3=np.r_[w3,0,0]
X3=np.fft.fft(x*w3)
w4=0.5+0.5*np.cos(2*np.pi*np.arange(-(N-2)/2+1,(N-2)/2)/(N-3))
w4=np.r_[w4,0,0,0]
X4=np.fft.fft(x*w4)
print 'DFTs (phases)'
print np.angle(X1[:5])
print np.angle(X2[:5])
print np.angle(X3[:5])
print np.angle(X4[:5])
print 'DFTs (magnitudes)'
print np.abs(X1[:5])
print np.abs(X2[:5])
print np.abs(X3[:5])
print np.abs(X4[:5])
print 'Angles at time of middle of window'
print np.angle(x[N/2])
print np.angle(x[N/2-1])
print np.angle(x[N/2-2])
print 'DFT with phase adjustment (phases)'
print np.angle((X1*np.exp(1j*2*np.pi*(N/2)/N)*n)[:5])
print np.angle((X2*np.exp(1j*2*np.pi*(N/2-1)/N)*n)[:5])
print np.angle((X3*np.exp(1j*2*np.pi*(N/2-1)/N)*n)[:5])
print np.angle((X4*np.exp(1j*2*np.pi*(N/2-2)/N)*n)[:5])
print 'DFT with phase adjustment (magnitudes)'
print np.abs((X1*np.exp(1j*2*np.pi*(N/2)/N)*n)[:5])
print np.abs((X2*np.exp(1j*2*np.pi*(N/2-1)/N)*n)[:5])
print np.abs((X3*np.exp(1j*2*np.pi*(N/2-1)/N)*n)[:5])
print np.abs((X4*np.exp(1j*2*np.pi*(N/2-2)/N)*n)[:5])
print 'DFT of rotated data (phases)'
X1r=x
X1r=np.roll(X1r,-N/2)
X1r=np.fft.fft(X1r*w1)
X2r=x
X2r=np.roll(X2r,-N/2+1)
X2r=np.fft.fft(X2r*w2)
X3r=x
X3r=np.roll(X3r,-N/2+1)
X3r=np.fft.fft(X3r*w3)
X4r=x
X4r=np.roll(X4r,-N/2+2)
X4r=np.fft.fft(X4r*w4)
print np.angle(X1r[:5])
print np.angle(X2r[:5])
print np.angle(X3r[:5])
print np.angle(X4r[:5])
print 'DFT of rotated data (magnitudes)'
print np.abs(X1r[:5])
print np.abs(X2r[:5])
print np.abs(X3r[:5])
print np.abs(X4r[:5])

