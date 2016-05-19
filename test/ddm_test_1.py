# A test computing the 2nd order argument parameters using DDM.
# WARNING: This program makes the computer lock up momentarily and report a
# MemoryError...
import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

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
    mi=mi[np.flatnonzero(np.greater(x[mi],th))]
    return mi

# Synthesize test signal
# For this file, the coefficients of the modulating polynomial are with respect
# to the sample number.

# Number of partials
P=1

# time in seconds
T=1.

# Sample rate
Fs=16000

# Length in samples
N=Fs*T

# Sample indices
n=np.arange(N)

# Range to look for matching maximum frequencies in
# looks in bins -de_K to de_K away from centre bin
de_K=2

# Polynomial coefficients

# real parts
a_r=np.array([
    # initial amplitude, small range of values
    np.random.uniform(1.e-1,1.,P),
    # 1st amplitude coefficient
    # picked so that amplitude does not go up more than by a factor of 2
    # and down by a factor of 100 from the beginning to the end of the signal
    0.1,#np.random.uniform(np.log(1.e-2),np.log(2.),P)/N,
    # The same but time is squared so divided by time squared
    0.1#np.random.uniform(np.log(1.e-2),np.log(2.),P)/(N**2.)
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
    0.5*np.random.uniform(-.1,.1,P)/(float(N))
    ])
a_i[2,:]*=a_i[1,:]
a_i=np.flipud(a_i)

x=np.exp(np.apply_along_axis(np.polyval,0,a_r,n)
        + complex('1j')*np.apply_along_axis(np.polyval,0,a_i,n))
xs=np.sum(x,1)

# Add noise
xs+=np.random.standard_normal(N)*np.sqrt(0.01)

# Size of DFT in STFT
N_stft=512

# Hop size for STFT
H_stft=128

# Length of the window
L_w=N_stft

# Window in STFT
#W_stft=np.hanning(N_stft)
W_0=0.5+0.5*np.cos(np.pi*2.*(np.arange(N_stft)-(L_w-1)/2.)/float(L_w))
dW_0=-np.pi/L_w*np.sin(2.*np.pi*(np.arange(N_stft)-(L_w-1)/2.)/float(L_w))

# frame up the signal to perform STFT
n_stft=np.add.outer(np.arange(0,N-N_stft,H_stft),np.arange(N_stft))
x_stft=xs[n_stft.astype('i')]
X_stft=x_stft*W_0
X_stft=np.fft.fft(X_stft)/np.sum(W_0)

# The centre frames
x_0=x_stft[1:-1,:]
# The frame one hop previous
x__H=x_stft[:-2,:]
# The frame one hop afterward
x_H=x_stft[2:,:]

# 1st polynomial coefficient
j1=complex('1j')
k_x=np.arange(N_stft)
X_0_dp_1_w=np.fft.fft(x_0*W_0)*np.exp(j1*2.*np.pi*(L_w-1)/2.*k_x)
X__H_dp_1_w=np.fft.fft(x__H*W_0)*np.exp(j1*2.*np.pi*(H_stft+(L_w-1)/2.)*k_x)
X_H_dp_1_w=np.fft.fft(x_H*W_0)*np.exp(j1*2.*np.pi*(-H_stft+(L_w-1)/2.)*k_x)

# 2nd polynomial coefficient
l_x=np.arange(-(L_w-1)/2,(L_w-1)/2+1)
l_x__H=l_x-H_stft
l_x_H=l_x+H_stft
X_0_dp_2_w=np.fft.fft(x_0*W_0*2.*l_x)*np.exp(j1*2.*np.pi*(L_w-1)/2.*k_x)
X__H_dp_2_w=np.fft.fft(x__H*W_0*2.*l_x__H)*np.exp(j1*2.*np.pi*(H_stft+(L_w-1)/2.)*k_x)
X_H_dp_2_w=np.fft.fft(x_H*W_0*2.*l_x_H)*np.exp(j1*2.*np.pi*(-H_stft+(L_w-1)/2.)*k_x)

# Derivative of window
X_0_wd=np.fft.fft(x_0*dW_0)#*np.exp(j1*2.*np.pi*(L_w-1)/2.*k_x)
X__H_wd=np.fft.fft(x__H*dW_0)#*np.exp(j1*2.*np.pi*(H_stft+(L_w-1)/2.)*k_x)
X_H_wd=np.fft.fft(x_H*dW_0)#*np.exp(j1*2.*np.pi*(-H_stft+(L_w-1)/2.)*k_x)

X0_dp1w_ma=[]
for row_ in np.abs(X_0_dp_1_w):
    zm_,zmi_=sm.lextrem(row_)
    X0_dp1w_ma.append(zmi_[lmaxw(zm_,12,6,1.e2)])

a_ddm=[]
for h in xrange(1,len(X0_dp1w_ma)-1):
    a_ddm.append([])
    for zmi_ in X0_dp1w_ma[h]:
        # don't consider maxima right on edge
        if (zmi_ < de_K) or (zmi_ > (N_stft-1-de_K)):
            continue
        pre_max=X__H_dp_1_w[h,zmi_-de_K:zmi_+de_K+1].argmax()+zmi_
        post_max=X_H_dp_1_w[h,zmi_-de_K:zmi_+de_K+1].argmax()+zmi_
        A=np.vstack((
            np.c_[X_0_dp_1_w[h,zmi_-1:zmi_+2],X_0_dp_2_w[h,zmi_-1:zmi_+2]],
            np.c_[X_H_dp_1_w[h,post_max-1:post_max+2],X_H_dp_2_w[h,post_max-1:post_max+2]],
            np.c_[X__H_dp_1_w[h,pre_max-1:pre_max+2],X__H_dp_2_w[h,pre_max-1:pre_max+2]]
        ))
        b=np.vstack((
            np.c_[X_0_dp_1_w[h,zmi_-1:zmi_+2]*-2.*np.pi*j1*np.arange(zmi_-1,zmi_+2)/L_w]
                +np.c_[X_0_wd[h,zmi_-1:zmi_+2]],
            np.c_[X_0_dp_1_w[h,post_max-1:post_max+2]*-2.*np.pi*j1*np.arange(post_max-1,post_max+2)/L_w]
                +np.c_[X_H_wd[h,post_max-1:post_max+2]],
            np.c_[X_0_dp_1_w[h,pre_max-1:pre_max+2]*-2.*np.pi*j1*np.arange(pre_max-1,pre_max+2)/L_w]
                +np.c_[X__H_wd[h,pre_max-1:pre_max+2]]
        ))
        a_ddm_=np.linalg.lstsq(A,-b)[0]
        a_ddm[-1].append(a_ddm_)

lmaxs=[]
maxis=[]
for row_ in np.abs(X_stft):
#    z_=np.nonzero((row_[1:] > row_[:-1]) & (row_[:-1] > row_[1:]))
    zm_,zmi_=sm.lextrem(row_)
    maxis.append(zmi_)
    lmaxs.append(lmaxw(zm_,12,6,1.e-1))
#    lmaxs.append(np.arange(len(zmi_)))

plt.figure(1)
plt.plot(n,np.abs(x))

plt.figure(2)
plt.imshow(np.abs(X_stft).T,origin='lower',interpolation='none')
#for l in xrange(len(lmaxs)):
#    plt.scatter([l for _ in xrange(len(lmaxs[l]))],maxis[l][lmaxs[l]],c='r')
for h in xrange(1,len(X0_dp1w_ma)-1):
    for a in a_ddm[h-1]:
        a1=np.imag(a[0])
        a2=np.imag(a[1])
        n0=-H_stft/2.
        n1=H_stft/2.
        w0=a1+n0*a2
        w1=a1+n1*a2
        k0=w0/(np.pi*2.)*L_w
        k1=w1/(np.pi*2.)*L_w
        plt.plot([n0/H_stft+h,n1/H_stft+h],[k0,k1],c='r')

plt.show()

