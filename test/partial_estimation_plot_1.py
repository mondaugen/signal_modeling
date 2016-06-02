import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys
import ptpath as pp
from cvxopt import solvers
from copy import deepcopy
import ptpath_test as ppt

# Read 1 channel float64 samples from standard in.
# Solve a linear program connecting a percentage of the minimum number of nodes
# in one of the frames' worth of paths.
# This one investigates the viability of splitting the STFT into blocks limited
# in both time and frequency

# Width of frequency block in bins
L_k = 100
# Hop size of frequency block in bins
H_k = 75

# This scales the minimum number of nodes in one frame of a block, and sets the
# number of paths to be discovered by the path searching algorithm
R_B=1
# Hop size in samples
H=1024
with open('tmp/ac_gtr_a3_op_sr16k.f64','r') as f:
    x=np.fromfile(f)
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
b_ddm_hz=70.
# Spacing between the start points of these bands (in Hz)
o_ddm_hz=35.
# Convert to bins
b_ddm=np.round(b_ddm_hz/Fs*N)
o_ddm=np.round(o_ddm_hz/Fs*N)
print 'b_ddm, o_ddm = %d %d' %(b_ddm,o_ddm)
# threshold of value seen as valid
th_ddm=0.01
# Highest bin to consider
M_ddm=N/2
# number of bins after last maximum to skip
i_ddm=4

def _cost_func(a,b):
    # The cost of predicting the frequency of b using the parameters of a
    return (((np.imag(a.value[1]) + np.imag(a.value[2])*H) -
            np.imag(b.value[1]))**2.) + 1.

while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,
        b_ddm,o_ddm,th_ddm,M_ddm,i_ddm))
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
a_e_db=-60

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
f_hops=np.arange(0.,M/2.-L_k,H_k)
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
    if (k_a_e-k_a < 2):
        continue
    S,F=tf_block(a,k_a,k_a_e-k_a,h_k,L_k,M)
    f_min=min([len(f_) for f_ in F])
    f_max=max([len(f_) for f_ in F])
    if (f_min < 1):
        continue
    N_nodes=len(S.keys())
    J_min=int(f_min*R_B)
    #J_min=int(f_min*float(H_k)/float(L_k)
    #        *float(N_nodes)/(float(N_nodes_tot)/float(N_f_hops)))
    if (J_min < 1):
        J_min = 1
    if (J_min > f_min):
        J_min=f_min
    print 'number of nodes: %d' % (len(S.keys()),)
    print 'bounds on num nodes in frame: min: %d, max: %d' % (f_min,f_max)
    C,C_cxn=ppt.shortest_paths_cost_lattice(S,F,J,_cost_func)
    q=ppt.shortest_paths_viterbi(C,C_cxn)
    #d=pp.g_f_2lp(S,F,J_min,_cost_func,{'calc_mean':0,'min_mean_dev':0})
    #solx_a.append(solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b'])['x'])
    K_a.append(k_a)
    K_a_e.append(k_a_e)
    S_a.append(S)
    F_a.append(F)

for s_a,f_a,x,k_a,k_a_e in zip(S_a,F_a,solx_a,K_a,K_a_e):
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

## extract paths
#paths=[]
#for k_a in xrange(len(a)-N_f+1):
#    paths.append([])
#    for s_a,f_a,x in zip(S_a[k_a],F_a[k_a],solx_a[k_a]):
#        paths[-1].append([])
#        # Follow the path according to the solution
#        N_nodes=len(s_a.keys())
#        for key in s_a.keys():
#            path=[]
#            key_=key
#            path.append(key_)
#            done=False
#            while not done:
#                k_=0
#                premature=False
#                while k_ < len(s_a[key_].out_nodes):
#                    nd=s_a[key_].out_nodes[k_]
#                    if x[key_+nd*N_nodes] > 0.5:
#                        key_=nd
#                        path.append(key_)
#                        premature=True
#                        break
#                    k_+=1
#                if not premature:
#                    # No connected node was found, end of path
#                    done=True
#            if (len(path) == N_f):
#                paths[-1][-1].append(path)
#
##print paths
#
#
#plt.figure(2)
## Compute the principal components of the paths
#for k_a in xrange(0,len(a)-N_f+1):
#    # compute number of paths.
#    npaths=sum([len(paths_) for paths_ in paths[k_a]])
#    # Allocate matrix of observations
##    X=np.ndarray((N_f*2,npaths))
#    X=np.ndarray((N_f,npaths))
#    col=0
#    for s_a,paths_ in zip(S_a[k_a],paths[k_a]):
#        for path in paths_:
#            for k in xrange(N_f):
#                a_=s_a[path[k]].value
##                X[k,col]=np.imag(a_[2])/np.imag(a_[1])
#                X[k,col]=np.real(a_[1])/np.imag(a_[1])
##                X[k*2,col]=np.real(a_[1])#/np.real(a_[0])
##                X[k*2+1,col]=np.imag(a_[2])/np.imag(a_[1])
#            col+=1
#    A=sm.pca_ne(X)
#    plt.scatter(k_a*np.ones(npaths),A[0,:],c='b')
#
## To prove it works, we plot the paths
#plt.figure(3)
#for k_a in xrange(0,len(a)-N_f+1,3):
#    for s_a,paths_ in zip(S_a[k_a],paths[k_a]):
#        for path in paths_:
#            for k in xrange(len(path)-1):
#                src=path[k]
#                dst=path[k+1]
#                n0=s_a[src].frame_num+k_a
#                n1=s_a[dst].frame_num+k_a
#                k0=np.imag(s_a[src].value[1])/(2.*np.pi)*M
#                k1=np.imag(s_a[dst].value[1])/(2.*np.pi)*M
#                plt.plot([n0,n1],[k0,k1],c='k')
#
#plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')

plt.show()
