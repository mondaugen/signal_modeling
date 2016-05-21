import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

H=128
N=4096
a_s=[]
a_s.append([
    complex('1.+2.j'),
    complex('1e-3+1e-3j'),
    complex('-3e-7+1e-4j'),
    complex('0-1e-8j')
    ])
a_s.append([
    complex('1.+2.j'),
   complex('1.5e-3+2.e-3j'),
    complex('-1e-6+1e-4j'),
    complex('0-1e-7j')
    ])
n=np.arange(N)
x=np.zeros(N).astype('complex_')
for a_ in a_s:
    x+=np.exp(a_[0]+a_[1]*n+a_[2]*n**2.+a_[3]*n**3.)
#x+=np.random.standard_normal(N)
y=np.zeros(N).astype('complex_')
plt.plot(n,np.real(x))
a=[]
h=0
M=512
m=np.arange(M)
w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
while ((h+M) <= N):
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,32,16,10.,M))
    for a_ in a[-1]:
        y[h:h+M]+=np.exp(a_[0]+a_[1]*m+a_[2]*m**2.)*w/2.
    h+=H
plt.plot(n,np.real(y))
plt.show()
