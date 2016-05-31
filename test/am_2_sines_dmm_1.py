# Extract amplitude modulation parameters using DDM of an amplitude modulated
# sinusoid.
# Try and fit a sinusoid to the amplitude envelope

import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

# A2
f1=110
# A#2
f2=116.54
Fs=16000
T=1.
N=T*Fs
n=np.arange(N)
x1=np.exp(1j*2.*np.pi*n/Fs*f1)
x2=np.cos(2.*np.pi*n/Fs*f2)
x=x1+x2+np.random.standard_normal(N)*np.sqrt(0.1)
t=n/Fs
plt.figure(1)
plt.plot(t,np.real(x))
N_w=512
H=128
n_w=np.arange(N_w)
# window coefficients
#wc=np.r_[0.358735,0.488305,0.141265,0.011695]
wc=np.r_[0.5,0.5]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./N_w*np.outer(np.arange(len(wc)),n_w)))
w=np.sum(w_,0)
dw_=((2.*np.pi/N_w*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
    *np.sin(np.pi*2./N_w*np.outer(np.arange(1,len(wc)),n_w)))
dw=np.sum(dw_,0)
a=[]
for h in np.arange(0,N-N_w,H):
#    a.append(sm.ddm_p2_3_3(x[h:h+2*H+N_w],H,w,dw,4))
    a.append(sm.ddm_p2_1_3(x[h:(h+N_w)],w,dw))
h_i=0
h=0
#x_a=np.zeros(N).astype('complex_')
x_a=np.zeros(N)
x_f=np.zeros(N)
for h in np.arange(0,N-N_w,H):
    a_=a[h_i]
#    x_a[h:h+N_w]+=np.exp(a_[0]+a_[1]*n_w+a_[2]*n_w**2.)*w/2.
    x_a[h:h+N_w]+=np.exp(np.real(a_[0]+a_[1]*n_w+a_[2]*n_w**2.))*w/2.
    x_f[h:h+H]=np.imag(a_[1]+2.*a_[2]*n_w[:H])/(np.pi*2.)*Fs
    h_i+=1
plt.plot(t,x_a,c='k')
plt.figure(2)
plt.plot(t,x_f,c='k')
plt.show()
