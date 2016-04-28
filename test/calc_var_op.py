# Calculate the variance for one of the variables of an optimization problem
import cvxopt as cvx
import numpy as np
from math import sqrt

N=10
cvx.setseed()
w=cvx.uniform(N)

c=cvx.matrix(0,(N+2,1),tc='d')
c[N]=1
A=cvx.matrix(0,(N+1,N+2),tc='d')
A[0,:N]=w.T/float(N)
A[0,N+1]=-1
A[1:,:N]=cvx.spdiag(cvx.matrix(1,(N,1),tc='d'))
b=cvx.matrix(1,(N+1,1),tc='d')
b[0]=0
G=cvx.matrix(0,(N+1,N+2),tc='d')
G[0,-2]=-1
G[1:,:-2]=-1*cvx.spdiag(w)/sqrt(float(N))
h=cvx.matrix(0,(N+1,1),tc='d')
dims={'l':0,'q':[N+1],'s':[]}
sol=cvx.solvers.conelp(c,G,h,dims,A,b)
print sol['x']
print 'Variance:'
print sol['x'][-2]**2.-sol['x'][-1]**2.
print 'True variance'
print np.var(np.array(list(w)))
