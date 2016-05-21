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
y=np.zeros(N)
plt.plot(n,np.real(x))
a=[]
h=0
M=512
m=np.arange(M)
w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
while ((h+M) <= N):
    a.append(sm.ddm_p2_1_3(x[h:(h+M)],w,dw))
    y[h:h+M]+=np.exp(a[-1][0]+a[-1][1]*m+a[-1][2]*m**2.)*w/2.
    h+=H
plt.plot(n,np.real(y))
plt.show()
