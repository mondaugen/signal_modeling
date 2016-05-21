import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys

H=512
#N=4096
x=np.fromfile(sys.stdin)
N=len(x)
n=np.arange(N)
y=np.zeros(N).astype('d')
#plt.plot(n,np.real(x))
a=[]
h=0
M=2048
m=np.arange(M)
w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,32,16,.01,M/2))
    for a_ in a[-1]:
        y[h:h+M]+=(np.exp(np.real(a_[0])+np.real(a_[1])*m+np.real(a_[2])*m**2.)
            *np.cos(np.imag(a_[0])+np.imag(a_[1])*m+np.imag(a_[2])*m**2.)*w)
    h+=H
#plt.plot(n,y)
#plt.show()
#y/=np.max(y)
y.tofile(sys.stdout)
