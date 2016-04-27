# Functions for finding best paths related to partial tracking
import numpy as np
import cvxopt as cvx
from math import sqrt

class LPNode:
    def __init__(self,value=0,in_nodes=[],out_nodes=[],frame_num=-1):
        self.value=value
        self.in_nodes=in_nodes
        self.out_nodes=out_nodes
        self.frame_num=frame_num

def LPNode_euclidean_dist(a,b):
    """Euclidean distance from a to b."""
    return sqrt(((a.value-b.value)**2.)+((a.frame_num-b.frame_num)**2.))

def g_f_2lp(S,F,J,cost_func=LPNode_euclidean_dist):
    """From a set of nodes representing a graph and a set of frames describing
    which nodes belong together, produce a set of matrices and vectors to solve
    a linear program for the solution of J best paths that do not overlap,
    unless they have to (because there are too few nodes in a frame).
        
        Input:
            S:  The set of nodes describing the graph.
            F:  A list containing tuples of node numbers representing the nodes
                in each frame.
            cost_func(S_i,S_j):
                A function accepting two nodes returning a value representing
                their distance, or some other cost of connecting them.

        Returns dictionary with fields:
            c:  a vector representing the distances between connected nodes.
                Unconnected nodes have the distance Inf.
            A:  the equality constraint matrix
            b:  the equality contraint vector
            G:  the inequality constraint matrix of the expression Gx<=h
            h:  the equality constraint matrix of the above expression
    """
    N_nodes=len(S)
    c=cvx.matrix(float('inf'),(N_nodes,N_nodes))
    cmax=float('-inf')
    for i in S.keys():
        for j in S[i].out_nodes:
            cst=cost_func(S[i],S[j])
            if cst > cmax:
                cmax=cst
            c[i,j]=cst
    rs,cs=c.size
    for r_ in xrange(rs):
        for c_ in xrange(cs):
            if (c[r_,c_]==float('inf')) or (c[r_,c_]==float('-inf')):
                c[r_,c_]=cmax+1

    c.size=(N_nodes*N_nodes,1)
    A_=list()
    b_=list()
    G_=list()
    h_=list()
    for f in F[:-1]:
        a1_=cvx.spmatrix([],[],[],(N_nodes,N_nodes),'d')
#        a2_=cvx.spmatrix([],[],[],(N_nodes,N_nodes),'d')
        for i in f:
            # restrict number of edges leaving S[i]
            a1_[i,:]=1
            a1_[i,i]=0
#            for j in S[i].out_nodes:
#                a1_[i,j]=1
            # restrict number of edges entering S[i]
#            for j in S[i].in_nodes:
#                a2_[j,i]=1
        A_.append(a1_)
#        A_.append(a2_)
        # restrict this to the number of paths
        b_.append(J)
#        b_.append(J)

    for k in xrange(1,len(F)-1):
        # For each node in the inner frames, the number of edges entering a node
        # must equal the number of edges exiting the same node
        for i in F[k]:
            a_=cvx.spmatrix([],[],[],(N_nodes,N_nodes),'d')
            for j in S[i].in_nodes:
                a_[j,i]=1
            for j in S[i].out_nodes:
                a_[i,j]=-1
            A_.append(a_)
            # this must sum to 0
            b_.append(0)

    for k in xrange(len(F)):
        if len(F[k]) >= J:
            # if number of nodes in frame is >= to the number of paths, each
            # node can have a maximum of 1 edge entering and 1 edge leaving
            for i in F[k]:
                if (k != (len(F)-1)):
                    # (in last frame, there can be no edges leaving)
                    g1_=cvx.spmatrix([],[],[],(N_nodes,N_nodes),'d')
                    # restrict number of edges leaving S[i]
                    for j in S[i].out_nodes:
                        g1_[i,j]=1
                    G_.append(g1_)
                    # restrict this to the number of paths
                    h_.append(1)
                if (k != 0):
                    # (in first frame, there can be no edges entering)
                    g2_=cvx.spmatrix([],[],[],(N_nodes,N_nodes),'d')
                    # restrict number of edges entering S[i]
                    for j in S[i].in_nodes:
                        g2_[j,i]=1
                    G_.append(g2_)
                    # restrict this to the number of paths
                    h_.append(1)

    # allocate equality constraint matrix
    A=cvx.spmatrix([],[],[],(len(A_),N_nodes*N_nodes),'d')
    # fill with values
    for i in xrange(len(A_)):
        # flatten matrix
        A_[i].size=(N_nodes*N_nodes,1)
        # store transpose (see LP definition)
        A[i,:]=A_[i].T
        # Unflatten
        A_[i].size=(N_nodes,N_nodes)
    # make equality contstraint vector
    b=cvx.matrix(b_,tc='d')

    # allocate inequality constraint matrix, including space for the constraints
    # on the variable (min 0, max 1)
    G=cvx.spmatrix([],[],[],(len(G_)+2*N_nodes*N_nodes,N_nodes*N_nodes),'d')
    for i in xrange(len(G_)):
        G_[i].size=(N_nodes*N_nodes,1)
        G[i,:]=G_[i].T
        G_[i].size=(N_nodes,N_nodes)
    for n in xrange(N_nodes*N_nodes):
        G[len(G_)+n,n]=1
        h_.append(1)
    for n in xrange(N_nodes*N_nodes):
        G[len(G_)+N_nodes*N_nodes+n,n]=-1
        h_.append(0)
    h=cvx.matrix(h_,tc='d')

    d=dict()
    d['c']=c
    d['A']=A
    d['A_']=A_
    d['b']=b
    d['b_']=b_
    d['G']=G
    d['G_']=G_
    d['h']=h
    d['h_']=h_
    return d
