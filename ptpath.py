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

def g_f_2lp(S,
            F,
            J,
            cost_func=LPNode_euclidean_dist,
            opt={'calc_mean':1,'min_mean_dev':1}):

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
            opt: 
                a dictionary with the following fields, set to 1 if the feature
                is desired
                calc_mean :
                    if set to 1, produces a set of constraints that will append
                    a vector of len(S) weights and len(F) means to the variable
                    x. The weights for each graph node are the sum of the
                    incoming edge activation coefficients and the outgoing edge
                    activation coefficients, divided by 2 (an average). The
                    means are the weighted average of node values (taken from
                    the node's "value" field) in a given
                    frame, divided by J (the node values are weighted by the
                    weights and then divided by J, giving the mean weighted mean
                    value)
                min_mean_dev :
                    (to use, calc_mean must also be set)
                    if set to 1, also adds variables that attempt to minimize
                    the deviation of the mean from 0 (by default). To make an
                    arbitrary value, set corresponding entries of the b
                    variable. You could also choose to weight this
                    proportionately to the costs determined by cost_func by
                    scaling entries of d['c'] accordingly.

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
        for i in f:
            # restrict number of edges leaving S[i]
            a1_[i,:]=1
            a1_[i,i]=0
        A_.append(a1_)
        # restrict this to the number of paths
        b_.append(J)

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

    A_mean=[]
#    b_mean=[]
    A_n_cols=N_nodes*N_nodes
    A_n_rows=len(A_)
    if (opt['calc_mean']):
        for k in xrange(len(F)):
            for i in F[k]:
                am1_=cvx.spmatrix([],[],[],(N_nodes*N_nodes+N_nodes+len(F),1),'d')
                if (k!=(len(F)-1)):
                    # For all frames but last frame, weight is simply sum of all
                    # edges exiting a node
                    # (indices are specified this way because the cost matrix is
                    # vectorized by stacking its columns)
                    for j in S[i].out_nodes:
                        am1_[N_nodes*j+i]=1
                if (k!=0):
                    # For all frames but first frame, weight is simply sum of all
                    # edges entering a node
                    for j in S[i].in_nodes:
                        am1_[N_nodes*i+j]=1
                if (k!=0) and (k!=(len(F)-1)):
                    # For all frames but first and last, weights are averaged
                    # between number of incoming and outgoing edges
                    am1_=0.5*am1_
                # Coefficient for variable that will contain this node's weight
                am1_[N_nodes*N_nodes+i]=-1
                A_mean.append(am1_)
                b_.append(0)
            # Make inequality for calculating mean from weights
            am2_=cvx.spmatrix([],[],[],(N_nodes*N_nodes+N_nodes+len(F),1),'d')
            for i in F[k]:
                am2_[N_nodes*N_nodes+i]=S[i].value/float(J)
            # Coefficient for variable that will contain the mean for this frame
            am2_[N_nodes*N_nodes+N_nodes+k]=-1
            A_mean.append(am2_)
            b_.append(0)
        A_n_cols+=N_nodes+len(F)
        A_n_rows+=len(A_mean)

    G_y=[]
    G_n_cols=len(G_)+2*N_nodes*N_nodes
    if (opt['min_mean_dev']):
        if (opt['calc_mean'] != 1):
            raise Exception('To use min_mean_dev, calc_mean must also be 1')
        for k in xrange(len(F)):
            gy1_=cvx.spmatrix([],[],[],(N_nodes*N_nodes+N_nodes+len(F)+len(F),1),'d')
            gy2_=cvx.spmatrix([],[],[],(N_nodes*N_nodes+N_nodes+len(F)+len(F),1),'d')
            gy3_=cvx.spmatrix([],[],[],(N_nodes*N_nodes+N_nodes+len(F)+len(F),1),'d')
            gy1_[N_nodes*N_nodes+N_nodes+k]=1
            gy1_[N_nodes*N_nodes+N_nodes+len(F)+k]=-1
            gy2_[N_nodes*N_nodes+N_nodes+k]=-1
            gy2_[N_nodes*N_nodes+N_nodes+len(F)+k]=-1
            gy3_[N_nodes*N_nodes+N_nodes+len(F)+k]=-1
            G_y.append(gy1_)
            G_y.append(gy2_)
            G_y.append(gy3_)
        G_n_cols+=len(G_y)
        A_n_cols+=len(F)

    # allocate equality constraint matrix
    A=cvx.spmatrix([],[],[],(A_n_rows,A_n_cols),'d')
    # fill with values
    A_idx=0
    for i in xrange(len(A_)):
        # flatten matrix
        A_[i].size=(N_nodes*N_nodes,1)
        # store transpose (see LP definition)
        A[i,:(N_nodes*N_nodes)]=A_[i].T
        # Unflatten
        A_[i].size=(N_nodes,N_nodes)
        A_idx+=1

    if (opt['calc_mean']):
        for i in xrange(len(A_mean)):
            A[A_idx,:(N_nodes*N_nodes+N_nodes+len(F))]=A_mean[i].T
            A_idx+=1

    # make equality contstraint vector
    b=cvx.matrix(b_,tc='d')

    # allocate inequality constraint matrix, including space for the constraints
    # on the variable (min 0, max 1)
    G=cvx.spmatrix([],[],[],(G_n_cols,A_n_cols),'d')
    G_idx=0
    for i in xrange(len(G_)):
        G_[i].size=(N_nodes*N_nodes,1)
        G[i,:(N_nodes*N_nodes)]=G_[i].T
        G_[i].size=(N_nodes,N_nodes)
        G_idx+=1
    for n in xrange(N_nodes*N_nodes):
        G[len(G_)+n,n]=1
        G_idx+=1
        h_.append(1)
    for n in xrange(N_nodes*N_nodes):
        G[len(G_)+N_nodes*N_nodes+n,n]=-1
        G_idx+=1
        h_.append(0)

    if (opt['min_mean_dev']):
        for i in xrange(len(G_y)):
            G[G_idx,:(N_nodes*N_nodes+N_nodes+len(F)+len(F))]=G_y[i].T
            G_idx+=1
            h_.append(0)

    h=cvx.matrix(h_,tc='d')
    
    c_final=cvx.matrix(0,(A_n_cols,1),tc='d')
    c_final[:(N_nodes*N_nodes)]=c

    d=dict()
    d['c']=c_final
    d['A']=A
    d['A_']=A_
    d['b']=b
    d['b_']=b_
    d['G']=G
    d['G_']=G_
    d['h']=h
    d['h_']=h_
    return d
