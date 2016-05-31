# Extract parameters with modified blackman window
import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

H=128
N=2048
a0=complex('1.+2.j')
a1=complex('1e-4+1e-3j')
a2=complex('-1e-6+1e-4j')
a3=complex('0.5e-9-1e-8j')
n=np.arange(N)
x=np.exp(a0+a1*n+a2*n**2.+a3*n**3.)
#x+=np.random.standard_normal(N)
y=np.zeros(N).astype('complex_')
plt.plot(n,np.real(x))
a=[]
h=0
M=512
m=np.arange(M)
w=0.5-0.5*np.cos(2.*np.pi*m/M)
dw=np.pi/M*np.sin(2.*np.pi*m/M)
wc=np.r_[0.358735,0.488305,0.141265,0.011695]
#wc=np.r_[0.5,0.5]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./M*np.outer(np.arange(len(wc)),m)))
w=np.sum(w_,0)
dw_=((2.*np.pi/M*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
    *np.sin(np.pi*2./M*np.outer(np.arange(1,len(wc)),m)))
dw=np.sum(dw_,0)
while ((h+M) <= N):
    a.append(sm.ddm_p2_1_3(x[h:(h+M)],w,dw))
    y[h:h+M]+=np.exp(a[-1][0]+a[-1][1]*m+a[-1][2]*m**2.)*w
    h+=H
plt.plot(n,np.real(y))
plt.show()
