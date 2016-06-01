# Here we estimate the amplitude modulations of one exponentially
# decaying sinusoid and see how well it works.
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

M=16000
N=1024
n=np.arange(N)
T60=2.
Fs=16000
f1=100.
H=256
m60=Fs*T60
m=np.arange(M)
a1_r=np.log(1.e-3)/m60
a1_i=2*np.pi*f1/Fs
a_r=a1_r*m
a_i=a1_i*np.ones(M)
a_r_e=np.zeros(M).astype('d')
a_i_e=np.zeros(M).astype('d')

x=np.exp((a1_r+1j*a1_i)*m)+np.random.standard_normal(M)*np.sqrt(.001)

wc=np.r_[0.358735,0.488305,0.141265,0.011695]
#wc=np.r_[0.5,0.5]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./N*np.outer(np.arange(len(wc)),n)))
w=np.sum(w_,0)
dw_=((2.*np.pi/N*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
    *np.sin(np.pi*2./N*np.outer(np.arange(1,len(wc)),n)))
dw=np.sum(dw_,0)
# Coefficient to make normalize overlapped windows
c_ol=(np.r_[w,np.zeros(H)]+np.r_[np.zeros(H),w]).max()

a=[]
x_=np.zeros(M).astype('complex_')
h=0
while h <= (M-N):
    a.append(sm.ddm_p2_1_3(x[h:h+M],w,dw))
    a_r_e[h:h+H]=np.real(a[-1][0]+a[-1][1]*n[:H])
    a_i_e[h:h+H]=np.imag(a[-1][1])
    x_[h:h+N]+=np.exp(a[-1][0]+a[-1][1]*n)*w/c_ol
    h+=H


mpl.rcParams['legend.fontsize'] = 10

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
ax.plot(m, a_i, a_r*20/np.log(10.), label='true', c='b')
ax.plot(m, a_i_e, a_r_e*20/np.log(10.), label='estimated', c='g')
ax.legend()

fig2=plt.figure(2)
ax=fig2.gca()
ax.plot(m,np.real(x),label='true')
ax.plot(m,np.real(x_),label='estimated')
ax.legend()

plt.show()
