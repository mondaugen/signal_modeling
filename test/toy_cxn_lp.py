from itertools import product
import matplotlib.pyplot as plt
from numpy import linalg, array
from cvxopt import matrix, spmatrix, solvers, sparse
import numpy as np
# Nodes in graph
N_nodes=15
nodes=[
        (1,7.4),
        (1,5.1),
        (1,3.8),
        (1,1.9),
        (2,5.7),
        (2,3.2),
        (3,7.6),
        (3,6.8),
        (3,5.1),
        (3,3.6),
        (3,2.1),
        (4,7.9),
        (4,6),
        (4,4.3),
        (4,3.5) ]
# Connected nodes
adjn=[  [(1,2,3,4),(5,6)],
        [(5,6),(7,8,9,10,11)],
        [(7,8,9,10,11),(12,13,14,15)]]
adjn=map(lambda z: list(product(z[0],z[1])),adjn)
C=matrix(100000,(N_nodes*N_nodes,1),'d')
for b in adjn:
    for a in b:
        i,j = a
        x=array(nodes[i-1])
        y=array(nodes[j-1])
        C[N_nodes*(i-1)+(j-1)]=linalg.norm(x-y)
#        plt.plot([x[0],y[0]],[x[1],y[1]],'k')
a_plt=array(nodes)
plt.scatter(a_plt[:,0],a_plt[:,1])
#plt.show()

# Put constraints on number of entering
E=dict()
e=dict()
# Number entering last frame
for a in [12,13,14,15]:
    a=a-1
    _m=np.zeros((N_nodes*N_nodes,1))
    _idx=N_nodes*a+np.arange(N_nodes)
    _m[_idx]=1
    E[a]=_m
    e[a]=1
# Number entering intermediate frames
for a in [(5,6),(7,8,9,10,11)]:
    a=np.array(a)-1
    _idx=(N_nodes*np.array(a,ndmin=2).T+np.arange(N_nodes))
    _m=np.zeros((N_nodes*N_nodes,1))
    _m[_idx]=1
    a=tuple(a)
    E[a]=_m
    e[a]=4
# Number of leaving
L=dict()
l=dict()
# Number leaving first frame
for a in [1,2,3,4]:
    a=a-1
    _idx=N_nodes*np.array(np.arange(N_nodes),ndmin=2).T+np.array(a)
    _m=np.zeros((N_nodes*N_nodes,1))
    _m[_idx]=1
    L[a]=_m
    l[a]=1
# number leaving intermediate frames
for a in [(5,6),(7,8,9,10,11)]:
    a=np.array(a)-1
    _idx=N_nodes*np.array(np.arange(N_nodes),ndmin=2).T+np.array(a)
    _m=np.zeros((N_nodes*N_nodes,1))
    _m[_idx]=1
    a=tuple(a)
    L[a]=_m
    l[a]=4
N_con=len(E)+len(L)
A=matrix(0,(N_nodes*N_nodes,N_con),'d')
b=matrix(0,(N_con,1),'d')
n=0
for k in E.keys():
    A[:,n]=matrix(E[k])
    b[n]=e[k]
    n+=1
for k in L.keys():
    A[:,n]=matrix(L[k])
    b[n]=l[k]
    n+=1
x_len=N_nodes*N_nodes
G=sparse([spmatrix(1,xrange(x_len),xrange(x_len)),-1*spmatrix(1,xrange(x_len),xrange(x_len))],'d')
h=matrix([matrix(1,(x_len,1)),matrix(0,(x_len,1))],tc='d')
print h.size
sol=solvers.lp(C,G,h,A=A.T,b=b)
print sol['status']
print sol['x']
s=matrix(sol['x'],(N_nodes,N_nodes))
Rw,Co=s.size
a_plt=array(nodes)
for r in xrange(Rw):
    for c in xrange(Co):
        if s[r,c] > 0.25:
            plt.plot(a_plt[(r,c),0],a_plt[(r,c),1],'k')
        
plt.show()
#plt.scatter(a_plt[:,0],a_plt[:,1])
#plt.show()
