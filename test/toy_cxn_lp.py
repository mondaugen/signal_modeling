from itertools import product
import matplotlib.pyplot as plt
from numpy import linalg, array
from cvxopt import matrix, spmatrix
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
C=matrix(float('inf'),(N_nodes*N_nodes,1))
for b in adjn:
    for a in b:
        i,j = a
        x=array(nodes[i-1])
        y=array(nodes[j-1])
        C[N_nodes*(i-1)+(j-1)]=linalg.norm(x-y)
        plt.plot([x[0],y[0]],[x[1],y[1]],'k')
#a_plt=array(nodes)
#plt.scatter(a_plt[:,0],a_plt[:,1])
#plt.show()

# Put constraints on number of entering
E=dict()
e=dict()
for a=[12,13,14,15]:
    _m=matrix(float('inf'),(N_nodes*N_nodes,1))
    _idx=N_nodes*a+np.arange(N_nodes)
    _m[_idx]=1
    E[a]=_m
    e[a]=1
