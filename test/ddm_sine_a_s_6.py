import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys

H=256
#N=4096
x=np.fromfile(sys.stdin)
N=len(x)
n=np.arange(N)
#plt.plot(n,np.real(x))
a=[]
h=0
M=1024
Y=np.zeros((M,(N-M)/H+1)).astype('complex_')
m=np.arange(M)
w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
k=0
w_s=np.sum(w)
while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,8,4,0.1,M/2))
    x_=x[h:(h+M)]
    Y[:,k]=np.fft.fft(x_*w)/w_s
    for a_ in a[-1]:
        w0=np.imag(a_[1])/(2.*np.pi)*M
        w1=(np.imag(a_[1])+2.*np.imag(a_[2])*(H-1))/(2.*np.pi)*M
        n0=k
        n1=k+1.*(H-1)/H
        plt.plot([n0,n1],[w0,w1],c='k')
    k+=1
    h+=H
plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')
plt.show()
#plt.plot(n,y)
#plt.show()
#y/=np.max(y)
#y.tofile(sys.stdout)
