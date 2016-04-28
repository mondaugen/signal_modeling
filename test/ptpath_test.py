# Define graph
from ptpath import LPNode, g_f_2lp, LPNode_euclidean_dist
import cvxopt as cvx
from cvxopt import solvers
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import random
import math
from math import sqrt, cos, sin, pi

S_test1={
        0: LPNode(0.12,[]     ,[2,3,4],0),
        1: LPNode(0.34,[]     ,[2,3,4],0),
        2: LPNode(4.56,[0,1]  ,[5,6,7],1),
        3: LPNode(2.78,[0,1]  ,[5,6,7],1),
        4: LPNode(0.9,[0,1]  ,[5,6,7],1),
        5: LPNode(6.23,[2,3,4],[8,9]  ,2),
        6: LPNode(3.45,[2,3,4],[8,9]  ,2),
        7: LPNode(1.67,[2,3,4],[8,9]  ,2),
        8: LPNode(5.89,[5,6,7],[]     ,3),
        9: LPNode(2.21,[5,6,7],[]     ,3)
        }
F_test1=[
        (0,1),
        (2,3,4),
        (5,6,7),
        (8,9)
        ]

def test_1():
    d=g_f_2lp(S_test1,F_test1,2)
    return d

def rand_graph_1(min_Fk,max_Fk,K):
    """Generate random graph of K frames with minimum min_Fk and maximum max_Fk
    vertices per frame. Each vertex has height between 0 and 1
    """
    grph=dict()
    n=0
    F=list()
    S=dict()
    for k in xrange(K):
        Fk=int(math.floor(random.random()*(float(max_Fk)-float(min_Fk)+1.)+float(min_Fk)))
        F.append(tuple(xrange(n,n+Fk)))
        n=n+Fk
    for k in xrange(K):
        f=F[k]
        r=sorted([random.random() for _ in xrange(len(f))])
        if (k > 0):
            in_nodes=list(F[k-1])
        else:
            in_nodes=[]
        if (k < (K-1)):
            out_nodes=list(F[k+1])
        else:
            out_nodes=[]
        for i,r_ in zip(f,r):
            S[i]=LPNode(r_,in_nodes,out_nodes,k)
    grph['S']=S
    grph['F']=F
    return grph

def sol_test_1(d):
    return solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b'])

def sol_test_mean_1(d,S,F):
    c=d['c']
    c[-len(F):]=1
    return solvers.lp(c,d['G'],d['h'],d['A'],d['b'])

def ranks_test_1(d):
    e=dict()
    e['A']=np.reshape(np.array(list(cvx.matrix(d['A']))),d['A'].size)
    e['G']=np.reshape(np.array(list(cvx.matrix(d['G']))),d['G'].size)
    e['B']=np.concatenate([e['A'],e['G']])
    e['A_rank']=linalg.matrix_rank(e['A'])
    e['G_rank']=linalg.matrix_rank(e['G'])
    e['B_rank']=linalg.matrix_rank(e['B'])
    return e

def plot_test(sol,S):
    for k in S.keys():
        plt.scatter(S[k].frame_num,S[k].value,c='k')
    x_=sol['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].frame_num,S[c_].frame_num,]),
                         np.array([S[r_].value,S[c_].value,]),'k')
    plt.show()

def plot_test_1(sol):
    for k in S_test1.keys():
        plt.scatter(S_test1[k].frame_num,S_test1[k].value,c='k')
    x_=sol['x']
    x_.size=(len(S_test1),len(S_test1))
    r,c=x_.size
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S_test1[r_].frame_num,S_test1[c_].frame_num,]),
                         np.array([S_test1[r_].value,S_test1[c_].value,]),'k')
    plt.show()

def plot_sin_1(sol,S,F,J):
    for k in S.keys():
        plt.scatter(S[k].frame_num,S[k].value,c='k')
    x_=sol['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].frame_num,S[c_].frame_num,]),
                         np.array([S[r_].value,S[c_].value,]),'k')
    x_tck=np.arange(0,(len(F)-1),0.1)
    for j_ in J:
        plt.plot(x_tck,j_[0]+j_[2]*np.cos(2.*pi*(j_[1]*x_tck+j_[3])),'g.')
    plt.show()

def plot_sin_compare_1(sol1,sol2,S,F,J):
    for k in S.keys():
        plt.scatter(S[k].frame_num,S[k].value,c='k')
        k0=S[k].frame_num-0.1
        k1=S[k].frame_num+0.1
        f0=S[k].value-0.1*S[k].fm/(2.*pi)
        f1=S[k].value+0.1*S[k].fm/(2.*pi)
        plt.plot([k0,k1],[f0,f1],c='Gray')
    x_=sol1['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].frame_num,S[c_].frame_num,]),
                         np.array([S[r_].value,S[c_].value,]),'k')
    x_=sol2['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].frame_num,S[c_].frame_num,]),
                         np.array([S[r_].value,S[c_].value,]),'b')
    x_tck=np.arange(0,(len(F)-1),0.1)
    for j_ in J:
        plt.plot(x_tck,j_[0]+j_[2]*np.cos(2.*pi*(j_[1]*x_tck+j_[3])),'g.')
    plt.show()

    
def plot_test_mean(sol,S,F):
    N_frames=len(F)
    for k in S.keys():
        plt.scatter(S[k].frame_num,S[k].value,c='k')
        plt.scatter(S[k].frame_num,sol['x'][(-N_frames+S[k].frame_num)],c='r',linewidths=[0])
    x_=sol['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].frame_num,S[c_].frame_num,]),
                         np.array([S[r_].value,S[c_].value,]),'k')
    plt.show()

def FMNode_dist(a,b):
    """a and b must be instances of LPNode class with a fm field (added in after
    instantiation for instance. This uses the frequency modulation as a measure
    of closeness of 2 nodes."""
    return 0.02*abs(b.fm-a.fm)+0.98*LPNode_euclidean_dist(a,b)
#    return LPNode_euclidean_dist(a,b)

def fm_graph_1(K,J,L=0):
    """
    Generate graph of K frames whose nodes represent len(J) underlying
    sinusoids. 
    J is a list of tuples containing the offset, the frequency, the amplitude,
    and the initial phase for each sinusoidal component. Frequency value is
    cycles / frame.  L is the number of spurious values added per frame
    (default 0)

    Returns:
        grph: Dictionary with fields
            F: list of node indices in each frame
            S: Dictionary containing the node instances
            
    """
    grph=dict()
    n=0
    F=list()
    S=dict()
    ph=[j_[3] for j_ in J]
    x=[0 for _ in J]
    dx=[0 for _ in J]
    for k in xrange(K):
        Fk=len(J)+L
        F.append(tuple(xrange(n,n+Fk)))
        n=n+Fk
    for k in xrange(K):
        for j in range(len(ph)):
            x[j]=J[j][0]+J[j][2]*cos(2.*pi*ph[j])
            dx[j]=-2.*pi*J[j][2]*sin(2.*pi*ph[j])
            ph[j]+=J[j][1]
        f=F[k]
        if (k > 0):
            in_nodes=list(F[k-1])
        else:
            in_nodes=[]
        if (k < (K-1)):
            out_nodes=list(F[k+1])
        else:
            out_nodes=[]
        for i,j in zip(f[:len(J)],range(len(ph))):
            S[i]=LPNode(x[j],in_nodes,out_nodes,k)
            S[i].fm=dx[j]
        rand_x=[random.random() for _ in xrange(L)]
        rand_dx=[random.gauss(0,sqrt(10.)) for _ in xrange(L)]
        for i,j in zip(f[len(J):],range(L)):
            S[i]=LPNode(rand_x[j],in_nodes,out_nodes,k)
            S[i].fm=rand_dx[j]
    grph['S']=S
    grph['F']=F
    return grph

def test_fm_1(grph,L_J):
    d=g_f_2lp(grph['S'],grph['F'],L_J,FMNode_dist,{'calc_mean':0,'min_mean_dev':0})
    return d

def test_fm_2(grph,L_J):
    d=g_f_2lp(grph['S'],grph['F'],L_J,opt={'calc_mean':0,'min_mean_dev':0})
    return d
