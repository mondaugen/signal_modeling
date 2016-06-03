import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys
import ptpath as pp
from cvxopt import solvers
from copy import deepcopy
import ptpath_test as ppt

# Read 1 channel float64 samples from standard in.
# Split into sub-bands in frequency
# Find first frame where peak exceeds threshold
# Find last frame where peak exceeds threshold
# Find paths that go from this frame to the end frame.

# In this implementation, the path search can skip over frames containing no
# nodes

# Width of frequency block in bins
L_k = 32
# Hop size of frequency block in bins
H_k = 24
# Bottom of first frequency block in bins
I_k = 5

# This scales the minimum number of nodes in one frame of a block, and sets the
# number of paths to be discovered by the path searching algorithm
R_B=1
# Hop size in samples
H=1024
with open('tmp/ac_gtr_a3_op_sr16k.f64','r') as f:
    x=np.fromfile(f)[:12000]
# Sample rate of file
Fs=16000
# Length of the file
N=len(x)
# sample indices of the file
n=np.arange(N)
# space to store the estimate ddm parameters
a=[]
# current hop
h=0
M=4096
Y=np.zeros((M,(N-M)/H+1)).astype('complex_')
m=np.arange(M)
# Windows for DDM
#w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
#dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
wc=np.r_[0.358735,0.488305,0.141265,0.011695]
#wc=np.r_[0.5,0.5]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./M*np.outer(np.arange(len(wc)),m)))
w=np.sum(w_,0)
dw_=((2.*np.pi/M*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
    *np.sin(np.pi*2./M*np.outer(np.arange(1,len(wc)),m)))
dw=np.sum(dw_,0)
# Coefficient to make normalized overlapped windows
c_ol=(np.r_[w,np.zeros(H)]+np.r_[np.zeros(H),w]).max()
k=0
w_s=np.sum(w)
k_S=0
plt.figure(1)

# Size of band over which local maximum is searched (in Hz)
b_ddm_hz=50.
# Spacing between the start points of these bands (in Hz)
o_ddm_hz=25.
# Convert to bins
b_ddm=np.round(b_ddm_hz/Fs*M)
o_ddm=np.round(o_ddm_hz/Fs*M)
print 'b_ddm, o_ddm = %d %d' %(b_ddm,o_ddm)
# threshold of value seen as valid
th_ddm=0.01
# Highest bin to consider
M_ddm=M/2
# number of bins after last maximum to skip
i_ddm=3
# ratio between local maximum and local minimum
ar_ddm=10

def _cost_func(a,b):
    # The cost of predicting the frequency of b using the parameters of a
#    return (((np.imag(a.value[1]) + np.imag(a.value[2])*H) -
#            np.imag(b.value[1]))**2.) + 1.
    return ((((np.imag(a.value[1]) + np.imag(a.value[2])*H) -
            np.imag(b.value[1]))**2.)*1. +
           ((np.real(a.value[0]) + np.real(a.value[1])*H +
               np.real(a.value[2])*H**2. - np.real(b.value[0]))**2.) + 1.)

while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,
        b_ddm,o_ddm,th_ddm,M_ddm,i_ddm,ar_ddm))
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

# For each sub-band in frequency, find the first frame where the peak power exceeds
# some threshold and the last frame where it goes below this threshold
a_s_db=-60
a_e_db=-110

# The solution for each block
solx_a=[]
# The nodes of each block
S_a=[]
# The node indices per frame of each block
F_a=[]
# The start frame indices of each block
K_a=[]
# The end frame indices of each block
K_a_e=[]
# Number of hops in frequency over the frames
f_hops=np.arange(I_k,M/2.-L_k,H_k)
N_f_hops=len(f_hops)

# iterate over bands
for h_k in f_hops:
    # Find the frame where the peak power exceeds a threshold
    k_a=0
    done=False
    while((not done) and (k_a < Y.shape[1])):
        print (20.*np.log10(np.abs(Y[h_k:h_k+L_k,k_a]))).max()
        if (20.*np.log10(np.abs(Y[:,k_a]))).max() > a_s_db:
            done=True
        else:
            k_a+=1
    if (k_a == Y.shape[1]):
        # No start frame was found, go to next frequency block
        continue
    print k_a
    k_a_e=k_a
    done=False
    while((not done) and (k_a_e < Y.shape[1])):
        if (20.*np.log10(np.abs(Y[h_k:h_k+L_k,k_a_e]))).max() < a_e_db:
            done=True
        else:
            k_a_e+=1
    if (k_a_e == Y.shape[1]):
        # If no end frame found, just uses the last frame
        k_a_e-=1
    print k_a_e
    if (k_a_e-k_a < 3):
        continue
    S,F=tf_block(a,k_a,k_a_e-k_a,h_k,L_k,M)
    f_min=min([len(f_) for f_ in F])
    f_max=max([len(f_) for f_ in F])
    if (f_min < 1):
        continue
    N_nodes=len(S.keys())
#    J_min=int(f_min*R_B)
#    #J_min=int(f_min*float(H_k)/float(L_k)
#    #        *float(N_nodes)/(float(N_nodes_tot)/float(N_f_hops)))
#    if (J_min < 1):
#        J_min = 1
#    if (J_min > f_min):
#        J_min=f_min
    J_min=1
    print 'number of nodes: %d' % (N_nodes),)
    print 'bounds on num nodes in frame: min: %d, max: %d' % (f_min,f_max)
    C,C_cxn=ppt.shortest_paths_cost_lattice(S,F,J_min,_cost_func)
    q=ppt.shortest_paths_viterbi(C,C_cxn)
    solx_a.append(ppt.sol_x_from_spv(S,F,q,C_cxn))
    #d=pp.g_f_2lp(S,F,J_min,_cost_func,{'calc_mean':0,'min_mean_dev':0})
    #solx_a.append(solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b'])['x'])
    K_a.append(k_a)
    K_a_e.append(k_a_e)
    S_a.append(S)
    F_a.append(F)

for s_a,f_a,x,k_a,k_a_e in zip(S_a,F_a,solx_a,K_a,K_a_e):
    N_nodes=len(s_a.keys())
    for n in xrange(len(x)):
        src=n%N_nodes
        dst=int(n/N_nodes)
        plt.scatter(s_a[src].frame_num+k_a,np.imag(s_a[src].value[1])/(2.*np.pi)*M,'b')
        if x[n] > 0.5:
            n0=s_a[src].frame_num+k_a
            n1=s_a[dst].frame_num+k_a
            k0=np.imag(s_a[src].value[1])/(2.*np.pi)*M
            k1=np.imag(s_a[dst].value[1])/(2.*np.pi)*M
            plt.plot([n0,n1],[k0,k1],c='k')

plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')

plt.show()
