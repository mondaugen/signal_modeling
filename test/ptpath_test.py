# Define graph
from ptpath import LPNode, g_f_2lp
import cvxopt as cvx
from cvxopt import solvers
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import random
import math

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
        print Fk
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
