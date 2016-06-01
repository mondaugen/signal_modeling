# Here we estimate the amplitude modulations of many exponentially
# decaying sinusoids and see how well it works.
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
# Fundamental frequency
# Be cautious of the window size, if the window is too short, the band width of
# 1 bin will be too great and harmonics too close to one another will distrupt
# the analysis.
f0=80.
# number of harmonics
K=10
k=np.arange(K)+1
f=f0*k
H=256
m60=Fs*T60
m=np.arange(M)
a1_r=np.log(1.e-3)/m60
a1_i=2*np.pi*f/Fs
a_r=np.log(k**-1.)[:,np.newaxis]+np.outer(np.ones(K),a1_r*m)
a_i=np.outer(a1_i,m)
a_r_e=[]
a_i_e=[]
b_ddm=4
o_ddm=2
th_ddm=2
M_ddm=N/2
i_ddm=1

#x=np.exp((a1_r+1j*a1_i)*m)+np.random.standard_normal(M)*np.sqrt(.001)
#x=np.exp((a1_r+1j*a1_i)*m)
x=np.sum(np.exp(a_r+1j*a_i),0)

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
    a.append(sm.ddm_p2_1_3_b(x[h:h+N],w,dw,
        b_ddm,o_ddm,th_ddm,M_ddm,i_ddm))
    a_r_e.append([])
    a_i_e.append([])
    for a_ in a[-1]:
        a_r_e[-1].append(np.real(a_[0]+a_[1]*n[:H]))
        a_i_e[-1].append(np.imag(a_[1])*np.ones(H))
        x_[h:h+N]+=np.exp(a_[0]+a_[1]*n)*w/c_ol
    h+=H

mpl.rcParams['legend.fontsize'] = 10

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
m_plt=np.outer(np.ones(K),m)
for k_ in xrange(K):
    ax.plot(m_plt[k_,1:], a_i[k_,1:]/m[1:], a_r[k_,1:]*20/np.log(10.), label='true', c='b')

h=0
for k in xrange(len(a_i_e)):
    for l in xrange(len(a_i_e[k])):
        ax.plot(np.arange(h,h+H).astype('d'),a_i_e[k][l],a_r_e[k][l]*20/np.log(10.),
            label='estimated',c='g')
    h+=H

#ax.plot(m, a_i_e, a_r_e*20/np.log(10.), label='estimated', c='g')
ax.legend()

fig2=plt.figure(2)
ax=fig2.gca()
ax.plot(m,np.real(x),label='true')
ax.plot(m,np.real(x_),label='estimated')
ax.legend()

plt.show()
