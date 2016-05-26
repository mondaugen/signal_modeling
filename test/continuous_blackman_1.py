import cvxopt
import numpy as np
import matplotlib.pyplot as plt
P=cvxopt.matrix(np.eye(4))
a_=np.array([.35875,.48829,.14128,.01168])
#a_=np.array([.40217,.49703,.09392,.00183])
q=-1.*cvxopt.matrix(a_)
A=cvxopt.matrix([[1,-1,1,-1],[1,1,1,1]],tc='d').T
b=cvxopt.matrix([0,1],tc='d')
a=cvxopt.solvers.qp(P,q,A=A,b=b)['x']
N=512
n=np.arange(N)
a=np.array(a)
print a.flatten()
w=np.sum(np.cos(np.outer(2.*np.pi*n/N,[0,1,2,3]))*np.power(-1,[0,1,2,3])*a.flatten(),
        axis=1)
#w=np.sum(np.cos(np.outer(2.*np.pi*n/N,[0,1,2,3]))*np.power(-1,[0,1,2,3])*a_.flatten(),
#        axis=1)
print w
print w[0],w[N-1]
W=np.fft.fft(w,16*N)/sum(w)
n_=np.arange(16*N)/16.
plt.figure(1)
plt.plot(n,w)
plt.figure(2)
plt.plot(n_,20*np.log10(np.abs(W)))
plt.show()
