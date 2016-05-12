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
import string
from itertools import combinations, product, permutations

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

def plot_test(sol,S,show=True,fignum=0):
    plt.figure(fignum)
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
    if (show):
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


    Returns a dictionary with following keys where
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
    # Principal components
    A=[]
    # Indices of values to keep
    Xdi=[]
    # Indices of values to discard
    Xki=[]
    # Classifications
    C=[]
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
        if (c.shape[1] == 8):
            # Load in principal components
            A.append(c[r,4])
            # load in keep and discard indices
            Xdi.append(c[r,5])
            Xki.append(c[r,6])
            # Load in classifications
            C.append(c[r,7])

        t_samp+=H
    rslt=dict()
    rslt['S']=S
    rslt['F']=F
    rslt['opt']=opt
    rslt['trues']=trues
    rslt['Xdi']=[x-1 for x in Xdi]
    rslt['Xki']=[x-1 for x in Xki]
    rslt['C']=C
    return rslt

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

def shortest_eu_path_viterbi(S,F):
    """
    Use the viterbi algorithm to find the best path through a trellis based on
    it having the minimum Euclidean distance.
    """
    F_lens=tuple([len(f) for f in F])
    T=len(F)
    S_=np.ndarray(1,dtype=string.join(['%dfloat64' for _ in
        xrange(T)],',')%F_lens)[0]
    de=np.ndarray(1,dtype=string.join(['%dfloat64' for _ in
        xrange(T)],',')%F_lens)[0]
    ps=np.ndarray(1,dtype=string.join(['%dint64' for _ in
        xrange(T)],',')%F_lens)[0]
    q=np.ndarray(T,'int64')
    de[0]=0
    ps[0]=0
    for t in xrange(T):
        for j in xrange(len(F[t])):
            S_[t][j]=S[F[t][j]].value

    for t in xrange(1,T):
        N=F_lens[t]
        for j in xrange(N):
            cst_=de[t-1]+np.power(S_[t][j]-S_[t-1],2)
            ps[t][j]=np.argmin(cst_)
            de[t][j]=cst_[ps[t][j]]

    q[T-1]=np.argmin(de[T-1])
    for t in (T-2-np.arange(T-1)):
        q[t]=ps[t+1][q[t+1]]
    return q

def plot_sepv(S,F,q):
    """
    Plot the shortest path returned by the viterbi algorithm.
    """
    T=len(F)
    ys=[]
    for t in xrange(T):
        for f in F[t]:
            plt.scatter(t,S[f].value,c='k')
        ys.append(S[F[t][q[t]]].value)
    plt.plot(list(xrange(T)),ys,c='g')

    plt.show()

def shortest_paths_cost_lattice(S,F,J,cost_func=LPNode_euclidean_dist):
    """
    Build the cost lattice for the Viterbi algorithm that find the best J
    non-overlapping paths through a lattice.
    The number of nodes in each frame must be >= to J.

    Returns (C,C_cxns)
        C: Are costs
        C_cxns: Are the connections assocated with the costs
    """
    T=len(F)
    # Cost of transitions
    C=[[] for _ in xrange(T-1)]
    # Cost of connections
    C_cxns=[[] for _ in xrange(T-1)]
    for t in xrange(T-1):
        it=product(combinations(xrange(len(F[t])),J),
                permutations(xrange(len(F[t+1])),J))
        for k in it:
            cst=0
            for i,j in zip(k[0],k[1]):
                cst+=cost_func(S[F[t+1][j]],S[F[t][i]])
            C[t].append(cst)
            C_cxns[t].append(tuple(k))
    return (C,C_cxns)

def shortest_paths_viterbi(C,C_cxn,big_cost=1000000):
    """
    Use viterbi algorithm to find shortest connected path through trellis C.
    """
    T=len(C)
    C_lens=tuple([len(c) for c in C])
    de=np.ndarray(1,dtype=string.join(['%dfloat64' for _ in
        xrange(T)],',')%C_lens)[0]
    ps=np.ndarray(1,dtype=string.join(['%dint64' for _ in
        xrange(T)],',')%C_lens)[0]
    q=np.ndarray(T,'int64')
    for j in xrange(C_lens[0]):
        de[0][j]=C[0][j]
    ps[0]=0

    for t in xrange(1,T):
        N=C_lens[t]
        for j in xrange(N):
            cst_=[0 for n_ in xrange(C_lens[t-1])]
            for i in xrange(C_lens[t-1]):
                if C_cxn[t-1][i][1] != C_cxn[t][j][0]:
                    # if the node indices don't match, give a prohibitively high
                    # cost
                    cst_[i]+=big_cost
                cst_[i]+=de[t-1][i]+C[t][j]
            ps[t][j]=np.argmin(cst_)
            de[t][j]=cst_[ps[t][j]]

    q[T-1]=np.argmin(de[T-1])
    for t in (T-2-np.arange(T-1)):
        q[t]=ps[t+1][q[t+1]]
    return q

def plot_spv(S,F,q,C_cxn,show=True,fignum=0):
    plt.figure(fignum)
    T=len(F)
    for t in xrange(T):
        for f in F[t]:
            plt.scatter(t,S[f].value,c='k')
    for t in xrange(T-1):
        for j,k in zip(C_cxn[t][q[t]][0],C_cxn[t][q[t]][1]):
            plt.plot([t,t+1],[S[F[t][j]].value,S[F[t+1][k]].value],'g')
    if (show):
        plt.show()

def plot_hsrpc_test(Z,D,show=True,fignum=0):
    """
    Using dictionary Z loaded from file, plot the different sources in different
    colours and the spurious sources in black.
    D is the dictionary created by reducing the nodes in Z['S'] (see
    ptpath_hsrp_cmp_1.py)
    """
    T=len(Z['Xki'])
    H=float(Z['opt'][0]['H'])
    plt.figure(fignum)
    for t in xrange(T):
        F_keep=np.array(Z['F'][t],dtype='i')
        F_keep=F_keep[Z['Xki'][t].astype('i')]
        for i in xrange(len(Z['C'][t])):
            w=float(D[int(F_keep[i])].value['w'])
            psi=float(D[int(F_keep[i])].value['psi'])
            w0=w-psi*H/2.
            w1=w+psi*H/2.
            plt.plot([t*H-H/2.,t*H+H/2.],
                    [w0,w1],
                    c=('#%06x' % (0xffffff/140*((int(Z['C'][t][i])+1)*30),)))
        F_dis=np.array(Z['F'][t],dtype='i')
        F_dis=F_dis[Z['Xdi'][t].astype('i')]
        for i in xrange(len(F_dis)):
            plt.scatter(t*H,
                    D[int(F_dis[i])].value['w'],
                    c='k')
    if(show):
        plt.show()

def make_pp_groups(Z,D):
    """
    From dictionary Z loded from file and dictionary created by reducing the
    nodes in Z['S'] (see plot_hsrpc_test), create a set S_groups containing
    nodes representing the data-points considered belonging to the same source
    in one frame, and also create a list of frames recording what nodes (i.e.,
    what keys in S_groups) are present in each frame.
    """
    T=len(Z['Xki'])
    F_groups=[[] for _ in xrange(T)]
    S_groups=dict()
    key=0
    for t in xrange(T):
        # The keys of the nodes in D that were kept (not considered spurious)
        F_keep=np.array(Z['F'][t],dtype='i')
        F_keep=F_keep[Z['Xki'][t].astype('i')]
        c_max=np.max(Z['C'][t])
        # iterate through the classifications
        for c in xrange(int(c_max)):
            S_t=list()
            for i in xrange(len(Z['C'][t])):
                if int(Z['C'][t][i]) == (c+1):
                    S_t.append(D[int(F_keep[i])].value)
            S_groups[key]=LPNode(S_t,[],[],t)
            F_groups[t].append(key)
            key+=1
        # Now associate the correct in and out nodes
    for t in xrange(T-1):
        # outnodes
        for f in F_groups[t]:
            S_groups[f].out_nodes=F_groups[t+1]
        # innodes
        for f in F_groups[t+1]:
            S_groups[f].in_nodes=F_groups[t]
    rslt=dict()
    rslt['S']=S_groups
    rslt['F']=F_groups
    return rslt

def LPNode_ppg_dist(a,b,opt):
    """
    Find the average error in predicting the values of b from the values of a.
    
    opt is a dictionary containing addtional information. In order to be able to
    pass this to the LP creation functions, you need to wrap it, like this:

    def my_func(a,b):
        opt={ ... } # stuff you want to pass
        return LPNode_ppg_dist(a,b,opt)

    ptpath.g_f_2_lp(...,cost_func=my_func)

    opt should contain the fields:
        H: the hop size in samples
        P(x,y): The probability of x xnd y being connected. Simply pass

                    def P(x,y):
                        return 1.

                if you want all to be considered equal. Otherwise pass a
                function like

                    def P(x,y):
                        exp(-(x['w']-y['w'])**2.)

                (it doesn't have to be a true probability, just a weighting)
    """
    t_a=a.frame_num
    t_b=b.frame_num
    t_de=(t_b-t_a)*opt['H']
    t_de=float(t_de)
    E_e=0
    for x in a.value:
        w_pred=x['w']+t_de*x['psi']
        for y in b.value:
            w_err=(y['w']-w_pred)**2.
            E_e+=w_err*opt['P'](x,y)
    return E_e

def lp_sol_extract_paths(sol,S,F):
    """
    Using a set of groups and frames, extract the connected paths.
    """
    T=len(F)
    N_F0=len(F[0])
    # Auxiliary variables will be tacked on the end, discard them
    x=sol['x'][:(len(S)*len(S))]
    x.size=(len(S),len(S))
    paths=[]
    for n in xrange(N_F0):
        # (we assume more than 1 frame)
        # Find connection for node S[F[0][n]]
        paths.append([])
        paths[-1].append(F[0][n])
        for t in xrange(1,T):
            K_Ft=len(F[t])
            for k in xrange(K_Ft):
                # If there was no path connection, the path still remains in the
                # paths list albeit with less nodes. We still check for
                # connections but there shouldn't be any.
                if x[paths[n][-1],F[t][k]] > 0.5:
                    paths[n].append(F[t][k])
    # Paths will now contain a list of paths represented by their keys in
    # S_groups (the set S containing the groups of nodes)
    return paths

def lp_sol_plot_paths(sol,S,F,show=True,fignum=0):
    """
    Plots the source paths for testing.
    """
    T=len(F)
    N_F0=len(F[0])
    x=sol['x'][:(len(S)*len(S))]
    x.size=(len(S),len(S))
    ro,co=x.size
    for k in S.keys():
        plt.scatter(S[k].value[0]['t_samp'],S[k].value[0]['w'],c='k')
    for r_ in xrange(ro):
        for c_ in xrange(co):
            if x[r_,c_] > 0.5:
                plt.plot([S[r_].value[0]['t_samp'],S[c_].value[0]['t_samp']],
                        [S[r_].value[0]['w'],S[c_].value[0]['w']],c='b')
    if(show):
        plt.show()


def pp_groups_plot_paths(S_groups,paths,opt,show=False,fignum=0):
    """
    Plot paths by looking up data in S_groups and plotting it in a color
    corresponding to the path.

    opt is a dictionary and should contain
        H: the hop size in samples
    """
    plt.figure(fignum)
    for c in xrange(len(paths)):
        t=0
        for p in paths[c]:
            for v in S_groups[p].value:
                t_=float(t*opt['H'])
                w_=float(v['w'])
                t0=-float(opt['H']/2.)
                t1=float(opt['H']/2.)
                w0=w_+t0*float(v['psi'])
                w1=w_+t1*float(v['psi'])
                plt.plot([t_+t0,t_+t1],[w0,w1],
                        c=('#%06x' % (0xffffff/140*((c+2)*30),)))
            t+=1

    if (show):
        plt.show()
