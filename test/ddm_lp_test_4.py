import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys
import ptpath as pp
from cvxopt import solvers
from copy import deepcopy

# Read 1 channel float64 samples from standard in.
# Solve a linear program connecting a percentage of the minimum number of nodes
# in one of the frames' worth of paths.
# This one investigates the viability of splitting the STFT into blocks limited
# in both time and frequency

# Argument 1 is the number of frames in a group for the LP, by default 4
N_f=6
if (len(sys.argv) >= 2):
    N_f=int(sys.argv[1])
# Width of frequency block in bins
L_k = 50.
# Hop size of frequency block
H_k = 37.5

R_B=0.65
H=256
x=np.fromfile(sys.stdin)
N=len(x)
n=np.arange(N)
a=[]
h=0
M=1024
Y=np.zeros((M,(N-M)/H+1)).astype('complex_')
m=np.arange(M)
# Windows for DDM
w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
k=0
w_s=np.sum(w)
k_S=0
plt.figure(1)

def _cost_func(a,b):
    # The cost of predicting the frequency of b using the parameters of a
    return (((np.imag(a.value[1]) + np.imag(a.value[2])*H) -
            np.imag(b.value[1]))**2.) + 1.

while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,6,3,0.01,M/2))
    x_=x[h:(h+M)]
    Y[:,k]=np.fft.fft(x_*w)/w_s
    k+=1
    h+=H

def tf_block(a,k_f,N_f,h_k,l_k,M):
    """
    a:
        The list of lists. Each leaf list contains vectors of peak parameters.
    k_f:
        The frame number.
    N_f:
        The number of frames to consider.
    h_k:
        The bottom frequency of the block.
    l_k:
        The length in the frequency dimension of the block.
    M:
        The total frequency dimension length.

    Returns

    S:
        The set of nodes in this block.
    F:  
        The node numbers in each frame of this block.
    """
    # Condense list
    a=a[k_f:k_f+N_f]
    a=[filter(lambda _a: ((np.imag(_a[1])/(2.*np.pi)*M >= h_k)
                and  (np.imag(_a[1])/(2.*np.pi)*M < (h_k+l_k))),
                a__) for a__ in a]
    S=dict()
    F=[]
    n_node=0
    for j in xrange(len(a)):
        F.append([])
        for k in xrange(len(a[j])):
            S[n_node]=pp.LPNode(a[j][k],[],[],j)
            F[-1].append(n_node)
            n_node+=1
    for j in xrange(1,len(F)-1):
        for k in xrange(len(F[j])):
            S[F[j][k]].in_nodes=F[j-1]
            S[F[j][k]].out_nodes=F[j+1]
    for k in xrange(len(F[0])):
        S[F[0][k]].out_nodes=F[1]
    for k in xrange(len(F[-1])):
        S[F[-1][k]].in_nodes=F[-2]
    return (S,F)

S_a=[]
F_a=[]
solx_a=[]

for k_a in xrange(len(a)-N_f+1):
    S_a.append([])
    F_a.append([])
    solx_a.append([])
    # Total number of nodes in these frames
    N_nodes_tot=sum([len(a_) for a_ in a[k_a:k_a+N_f]])
    # Number of hops in frequency over the frames
    f_hops=np.arange(0.,M/2.-L_k,H_k)
    N_f_hops=len(f_hops)
    for h_k in f_hops:
        S,F=tf_block(a,k_a,N_f,h_k,L_k,M)
        f_min=min([len(f_) for f_ in F])
        if (f_min < 1):
            continue
        S_a[-1].append(S)
        F_a[-1].append(F)
        N_nodes=len(S.keys())
        J_min=int(f_min*R_B)
        #J_min=int(f_min*float(H_k)/float(L_k)
        #        *float(N_nodes)/(float(N_nodes_tot)/float(N_f_hops)))
        if (J_min < 1):
            J_min = 1
        if (J_min > f_min):
            J_min=f_min
        d=pp.g_f_2lp(S,F,J_min,_cost_func,{'calc_mean':0,'min_mean_dev':0})
        solx_a[-1].append(solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b'])['x'])

for k_a in xrange(len(a)-N_f+1):
    for s_a,f_a,x in zip(S_a[k_a],F_a[k_a],solx_a[k_a]):
        N_nodes=len(s_a.keys())
        for n in xrange(len(x)):
            if x[n] > 0.5:
                src=n%N_nodes
                dst=int(n/N_nodes)
                n0=s_a[src].frame_num+k_a
                n1=s_a[dst].frame_num+k_a
                k0=np.imag(s_a[src].value[1])/(2.*np.pi)*M
                k1=np.imag(s_a[dst].value[1])/(2.*np.pi)*M
                plt.plot([n0,n1],[k0,k1],c='k')

plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')

# extract paths
paths=[]
for k_a in xrange(len(a)-N_f+1):
    paths.append([])
    for s_a,f_a,x in zip(S_a[k_a],F_a[k_a],solx_a[k_a]):
        paths[-1].append([])
        # Follow the path according to the solution
        N_nodes=len(s_a.keys())
        for key in s_a.keys():
            path=[]
            key_=key
            path.append(key_)
            done=False
            while not done:
                k_=0
                premature=False
                while k_ < len(s_a[key_].out_nodes):
                    nd=s_a[key_].out_nodes[k_]
                    if x[key_+nd*N_nodes] > 0.5:
                        key_=nd
                        path.append(key_)
                        premature=True
                        break
                    k_+=1
                if not premature:
                    # No connected node was found, end of path
                    done=True
            paths[-1][-1].append(path)

#print paths

# To prove it works, we plot the paths
plt.figure(2)
for k_a in xrange(0,len(a)-N_f+1,3):
    for s_a,paths_ in zip(S_a[k_a],paths[k_a]):
        for path in paths_:
            for k in xrange(len(path)-1):
                src=path[k]
                dst=path[k+1]
                n0=s_a[src].frame_num+k_a
                n1=s_a[dst].frame_num+k_a
                k0=np.imag(s_a[src].value[1])/(2.*np.pi)*M
                k1=np.imag(s_a[dst].value[1])/(2.*np.pi)*M
                plt.plot([n0,n1],[k0,k1],c='k')

plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')
plt.show()
