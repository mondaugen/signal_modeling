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
import scipy.io as sio

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
    d=g_f_2lp(S_test1,F_test1,2,opt={'calc_mean':0,'min_mean_dev':0})
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

def load_grph_harm_sines_rp(fpath):
    """
    Load data from fpath representing nodes in the graph but with more
    information to allow for more sophisticated cost functions. In particular,
    the data in the first column of the cell array are the grouping parameters
    (frequency modulation / frequency and amplitude modulation) and the data in
    the second column are matrices whose rows are the values in this order:
        [w,psi,abs(X),arg(X),mu]
    (see rm.m for a description)
    The data in the third column of the array are structures indicating the
    synthesis parameters. From here we are interested in the hop size so we know
    how separate different nodes are in time.
    This information is parsed and a graph is created which can then be passed
    to g_f_2lp along with a suitable cost calculation function to create a
    linear program to resolve trajectories in the graph. Note: before passing S
    to g_f_2lp you must reduce it to a dictionary indicating the node numbers
    rather than a list of lists by doing
    
        S_=reduce(lambda x,y: x+y,S)
        F_=reduce(lambda x,y: x+y,F)
        D=dict()
        for s,f in zip(S_,F_):
            D[f]=s


    Returns the tuple (S,F) where
        S: is the set of nodes describing the graph
        F: a list of tuples of node numbers describing the nodes numbers in each
           frame.
        opt: the contents of the 3rd column of the cell array. It can be used to
             synthesize the true parameters.
        trues: list of dictionaries containing the true values (not associated
               with any nodes though)
    """
    mat_contents=sio.loadmat(fpath)
    c=mat_contents['out']
    last_nodes=[]
    S=[]
    F=[]
    rows=c.shape[0]
    node_n=0
    H=float(c[0,2]['H'])
    opt=c[:,2]
    # time in samples
    t_samp=0
    trues=[]
    for r in xrange(rows):
        new_nodes=[]
        # combine first and second column matricies
        _rows=c[r,0].shape[0]
        S.append(list())
        for _r in xrange(_rows):
            M=np.hstack((c[r,0][_r,:],c[r,1][_r,:]))
            d={
                'psi_w' : M[0],
                'mu'    : M[1],
                'w'     : M[2],
                'psi'   : M[3],
                'a'     : M[4],
                'phi'   : M[5],
                't_samp': t_samp
            }
            s=LPNode(d,last_nodes,[],r)
            new_nodes.append(node_n)
            node_n+=1
            S[r].append(s)
        F.append(new_nodes)
        if (r>0):
            for _n in xrange(len(S[r-1])):
                S[r-1][_n].out_nodes=new_nodes
        last_nodes=new_nodes
        # now load in "true" values
        _rows=c[r,3].shape[0]
        for _r in xrange(_rows):
            M=c[r,3][_r,:]
            d={
                'mu'    : M[4],
                'w'     : M[0],
                'psi'   : M[1],
                'a'     : M[2],
                'phi'   : M[3],
                't_samp': t_samp
            }
            trues.append(d)
        t_samp+=H
    return (S,F,opt,trues)

def LPNode_rp_dist(a,b):
    """
    Calculates distance from a to b by extrapolating the w parameter of a using
    psi and comparing it to b's w parameter.
    """
    t1=b.value['t_samp'] - a.value['t_samp']
    if (t1 > 1000):
        t1=1000
    if (t1 < -1000):
        t1=-1000
    wb_=a.value['w']+a.value['psi']*t1
    # add 1 because we were have singularity issues with small values and all
    # values are positive so it shouldn't make a difference
    return (b.value['w'] - wb_) ** 2. + 1

def LPNode_rp_eu_dist(a,b):
    """
    Calculates distance from a to b as the euclidean distance from a.value['w']
    tp b.value['w'].
    """
    return (b.value['w'] - a.value['w']) ** 2.

def plot_lp_hsrp(sol,S,trues):
    for k in S.keys():
        plt.scatter(S[k].value['t_samp'],S[k].value['w'],c='k')

    x_=sol['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].value['t_samp'],S[c_].value['t_samp']]),
                         np.array([S[r_].value['w'],S[c_].value['w']]),'k')

    for t in trues:
        plt.scatter(t['t_samp'],t['w'],c='g')
    
    plt.show()

def plot_lp_hsrp_cmp(sol1,sol2,S,trues):
    plt.figure(1)
    xmin=0
    xmax=0
    ymin=0
    ymax=0
    for k in S.keys():
        plt.scatter(S[k].value['t_samp'],S[k].value['w'],c='k')
        if S[k].value['t_samp'] < xmin:
            xmin=S[k].value['t_samp']
        if S[k].value['t_samp'] > xmax:
            xmax=S[k].value['t_samp']
        if S[k].value['w'] < ymin:
            ymin=S[k].value['w']
        if S[k].value['w'] > ymax:
            ymax=S[k].value['w']
    xdiff=xmax-xmin
    ydiff=ymax-ymin
    plt.figure(2)
    for k in S.keys():
        plt.scatter(S[k].value['t_samp'],S[k].value['w'],c='k')

    x_=sol1['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    plt.figure(1)
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].value['t_samp'],S[c_].value['t_samp']]),
                         np.array([S[r_].value['w'],S[c_].value['w']]),'k')

    x_=sol2['x'][:(len(S)*len(S))]
    x_.size=(len(S),len(S))
    r,c=x_.size
    plt.figure(2)
    for r_ in xrange(r):
        for c_ in xrange(c):
            if x_[r_,c_] > 0.5:
                plt.plot(np.array([S[r_].value['t_samp'],S[c_].value['t_samp']]),
                         np.array([S[r_].value['w'],S[c_].value['w']]),'k')

    plt.figure(1)
    for t in trues:
        plt.scatter(t['t_samp'],t['w'],c='g',edgecolor='face')
    plt.title('Partial trajectories using frequency prediction error')
    plt.xlabel('Sample Number')
    plt.ylabel('Frequency in radians per sample')
    plt.figure(2)
    plt.title('Partial trajectories using Euclidean distance between frequencies')
    plt.xlabel('Sample Number')
    plt.ylabel('Frequency in radians per sample')
    for t in trues:
        plt.scatter(t['t_samp'],t['w'],c='g',edgecolor='face')

    plt.figure(1,figsize=(7.,4.5))
    plt.xlim((xmin-0.025*xdiff,xmax+0.025*xdiff))
    plt.ylim((ymin-0.025*ydiff,ymax+0.025*ydiff))
    plt.savefig('/tmp/lp_ptrack_fpredict.eps',bbox_inches='tight')
    plt.figure(2,figsize=(7.,4.5))
    plt.xlim((xmin-0.025*xdiff,xmax+0.025*xdiff))
    plt.ylim((ymin-0.025*ydiff,ymax+0.025*ydiff))
    plt.savefig('/tmp/lp_ptrack_eudist.eps',bbox_inches='tight')
    
    plt.show()
