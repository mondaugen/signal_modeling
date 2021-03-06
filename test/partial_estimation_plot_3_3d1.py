import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys
import ptpath as pp
from cvxopt import solvers
from copy import deepcopy
import ptpath_test as ppt
import pickle

# Read 1 channel float64 samples from standard in.
# Split into sub-bands in frequency
# Find first frame where peak exceeds threshold
# Find first frame from end where peak exceeds threshold
# Find paths that go from this frame to the end frame.

# In this implementation, the path search can skip over frames containing no
# nodes


# Width of frequency block in bins, this is the width of the window used when
# searching for peaks
L_k = 8 
# The width scalar. This is the percentage of the peak search window used when
# restricting what peaks are used in a partial path search. 1 means half width
# of window on either side of peak, 0.5 means quarter width of window on either
# side of peak, etc.
L_k_C=0.5
# Hop size of frequency block in bins
H_k = 4
# Bottom of first frequency block in bins
I_k = 50
# The number of bins to ignore after a maximum
i_k = 2
# The threshold of a true peak
th_k=1e-5
# The ratio of peak to its local minima (in dB)
ar_k=10

# A scalar adjusting the mean cost under which partial paths are set. A number
# greater than 1 will let more paths under the threshold, a number less than 1
# will let less
mc_C=1.2
# the path length at which the path thresholding function reaches 1.01 times the
# value mc*mc_C
mc_x0=50

# This scales the minimum number of nodes in one frame of a block, and sets the
# number of paths to be discovered by the path searching algorithm
R_B=1
# Hop size in samples
H=512

# Where to save partial trajectories
fout='tmp/xylo_fs4_sr16k.pardat'
with open('tmp/xylo_fs4_sr16k.f64','r') as f:
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
b_ddm_hz=15.
# Spacing between the start points of these bands (in Hz)
o_ddm_hz=7.5
# Convert to bins
b_ddm=np.round(float(b_ddm_hz)/Fs*M)
o_ddm=np.round(float(o_ddm_hz)/Fs*M)
print 'b_ddm, o_ddm = %d %d' %(b_ddm,o_ddm)
# threshold of value seen as valid
th_ddm=0.0001
# Highest bin to consider
M_ddm=M/2
# number of bins after last maximum to skip
i_ddm=2
# ratio between local maximum and local minimum
ar_ddm=3

def _cost_func(a,b):
    # The cost of predicting the frequency of b using the parameters of a
#    return (((np.imag(a.value[1]) + np.imag(a.value[2])*H) -
#            np.imag(b.value[1]))**2.) + 1.
#    return ((((np.imag(a.value[1]) + np.imag(a.value[2])*H) -
#            np.imag(b.value[1]))**2.)*1. +
#           ((np.real(a.value[0]) + np.real(a.value[1])*H +
#               np.real(a.value[2])*H**2. - np.real(b.value[0]))**2.) + 1.)
    return (np.imag(a.value[1]) - np.imag(b.value[1]))**2. + 1.

while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,
        b_ddm,o_ddm,th_ddm,M_ddm,i_ddm,ar_ddm))
    x_=x[h:(h+M)]
    Y[:,k]=np.fft.fft(x_*w)/w_s
    k+=1
    h+=H

def tf_block(a,k_f,N_f,h_k,l_k,M,a_rem):
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
    a_rem:
        A list of lists the same size as a. If element at same indices of a is
        True, then that data point has not yet been used.

    Returns

    S:
        The set of nodes in this block.
    F:  
        The node numbers in each frame of this block.
    """
    # Condense list
    a=a[k_f:k_f+N_f]
    a_rem=a_rem[k_f:k_f+N_f]
#    a=[filter(lambda _a: ((np.imag(_a[1])/(2.*np.pi)*M >= h_k)
#                and  (np.imag(_a[1])/(2.*np.pi)*M < (h_k+l_k))),
#                a__) for a__ in a]
    def _in_freq_bounds(_a):
        return ((np.imag(_a[1])/(2.*np.pi)*M >= h_k)
                and  (np.imag(_a[1])/(2.*np.pi)*M < (h_k+l_k)))

    a_=[]
    for j_ in xrange(len(a)):
        a_.append([])
        for k_ in xrange(len(a[j_])):
            if (_in_freq_bounds(a[j_][k_]) and a_rem[j_][k_]):
                a_[-1].append(a[j_][k_])
    a=a_

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
a_s_db=-100
a_e_db=-100

# The solution for each block
solx_a=[]
# The cost of the solution for each block
solc_a=[]
# A list of lists of indices of indices in each frame that are on the path
q_a=[]
# The nodes of each block
S_a=[]
# The node indices per frame of each block
F_a=[]
# The start frame indices of each block
K_a=[]
# The end frame indices of each block
K_a_e=[]
# The nodes that have not been classified into a path
a_rem=[[True for a__ in a_] for a_ in a]
## Number of hops in frequency over the frames
#f_hops=np.arange(I_k,M/2.-L_k,H_k)
#N_f_hops=len(f_hops)

# Find local maxima in initial spectrum to set f_hops
I_y_s=0
I_y_e=3
Y_ini=np.sum(np.abs(Y[:,I_y_s:I_y_e]),1)*(1./(I_y_e-I_y_s))
print Y_ini
f_hops=sm.lextrem_win(Y_ini,L_k,H_k,i_k,th_k,ar_k,M/2,I_k)


# iterate over bands
print f_hops
for h_k in f_hops:
    # Find the frame where the peak power exceeds a threshold
    k_a=0
    done=False
    # range of bins to consider
    r_h_k_l=int(h_k-np.round(L_k/2.*L_k_C))
    # if less than 0, set 0
    if (r_h_k_l < 0):
        r_h_k_l=0
    r_h_k_r=int(h_k+np.round(L_k/2.*L_k_C)+1)
    # if greater than max bin considered, set max bin considered
    if (r_h_k_r > M/2):
        r_h_k_r=M/2
    r_h_k=np.arange(r_h_k_l,r_h_k_r)
    while((not done) and (k_a < Y.shape[1])):
        print (20.*np.log10(np.abs(Y[r_h_k,k_a]))).max()
        if (20.*np.log10(np.abs(Y[:,k_a]))).max() > a_s_db:
            done=True
        else:
            k_a+=1
    if (k_a == Y.shape[1]):
        # No start frame was found, go to next frequency block
        continue
    print k_a
    k_a_e=Y.shape[1]-1
    done=False
    # Find the first frame from the end where peak power exceeds a threshold
    while((not done) and (k_a_e > k_a)):
        if (20.*np.log10(np.abs(Y[r_h_k,k_a_e]))).max() > a_e_db:
            done=True
        else:
            k_a_e-=1
    print k_a_e
    if ((k_a_e - k_a) < 2):
        continue
    print r_h_k_l
    print r_h_k_r
    S,F=tf_block(a,k_a,k_a_e-k_a,r_h_k_l,r_h_k_r-r_h_k_l,M,a_rem)
    print 'number of nodes: %d' % (len(S.keys()),)
    F=filter(len,F)
    if (len(F) < 3):
        continue
    f_min=min([len(f_) for f_ in F])
    f_max=max([len(f_) for f_ in F])
#    if (f_min < 1):
#        continue
    N_nodes=len(S.keys())
    J_min=int(f_min*R_B)
#    #J_min=int(f_min*float(H_k)/float(L_k)
#    #        *float(N_nodes)/(float(N_nodes_tot)/float(N_f_hops)))
    if (J_min < 1):
        J_min = 1
    if (J_min > f_min):
        J_min=f_min
    J_min=1
    print 'bounds on num nodes in frame: min: %d, max: %d' % (f_min,f_max)
    C,C_cxn=ppt.shortest_paths_cost_lattice(S,F,J_min,_cost_func)
    (q,c_tot)=ppt.shortest_paths_viterbi(C,C_cxn)
    solx_a.append(ppt.sol_x_from_spv(S,F,q,C_cxn))
    solc_a.append(c_tot)
    q_a.append(q)
    # indices of solution
    ## Indicate nodes used in this solution in order to not use them again
    #idx_=ppt.sol_S_from_spv(S,F,q,C_cxn)
    #for i_ in idx_:
    #    for k_ in xrange(len(a)):
    #        for l_ in xrange(len(a[k_])):
    #            if (a[k_][l_] is S[i_].value):
    #                a_rem[k_][l_]=False

    #d=pp.g_f_2lp(S,F,J_min,_cost_func,{'calc_mean':0,'min_mean_dev':0})
    #solx_a.append(solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b'])['x'])
    K_a.append(k_a)
    K_a_e.append(k_a_e)
    S_a.append(S)
    F_a.append(F)

mc=0.
c_q=[]
for q_,c_ in zip(q_a,solc_a):
    c_q.append(c_/len(q_))
    mc+=c_/len(q_)
mc/=len(solc_a)

mc_adj=(mc-1.)*mc_C+1.
mc_y_max=max(c_q)-1.
mc_b=np.log(mc_y_max)
mc_a=(np.log(0.01*mc_y_max)-mc_b)/float(mc_x0)

print mc
print mc_adj

f_avg=[]
da_avg=[]

# Information associated with each partial
p_info=[]
for s_a,f_a,x,k_a,k_a_e,cq_,q_ in zip(S_a,F_a,solx_a,K_a,K_a_e,c_q,q_a):
    # only plot if cost/length ratio <= the mean cost/length ratio
    # If a kept path, compute its average frequency
    _f_avg=0.
    # Variables to store amplitude of first and last data point
    _da_avg=0.
    _p_info=[]
    if(cq_ <= (np.exp(mc_a*len(q_)+mc_b)+mc_adj)):
        N_nodes=len(s_a.keys())
        for n in xrange(len(x)):
            src=n%N_nodes
            dst=int(n/N_nodes)
            if x[n] > 0.5:
                n0=s_a[src].frame_num+k_a
                n1=s_a[dst].frame_num+k_a
                k0=np.imag(s_a[src].value[1])/(2.*np.pi)*M
                k1=np.imag(s_a[dst].value[1])/(2.*np.pi)*M
                plt.plot([n0,n1],[k0,k1],c='k')
                _f_avg+=k0
                _da_avg+=np.real(s_a[src].value[1])
                _p_info.append([n0*H,s_a[src].value])
    else:
        # otherwise set its average frequency to -1 (ignored).
        _f_avg=-1.*len(q_)
    f_avg.append(float(_f_avg/len(q_)))
    da_avg.append(_da_avg/len(q_))
    p_info.append(_p_info)

# Save partial information
with open(fout,'w') as f:
    pickle.dump(p_info,f)

#for h_k in f_hops:
#    plt.plot([0,Y.shape[1]],[h_k,h_k],'g')
#    plt.plot([0,Y.shape[1]],[h_k+L_k,h_k+L_k],'r')
plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')

plt.plot(np.zeros(len(f_hops)),f_hops,'bo')

# plot path cost versus length
plt.figure(2)
print len(q_a),len(solc_a)
for q_,c_ in zip(q_a,solc_a):
    plt.plot(len(q_),c_/len(q_),'b.')

#plt.plot([0,max([len(q_) for q_ in q_a])],[mc_adj,mc_adj])
mc_plt_x=np.arange(0,max([len(q_) for q_ in q_a]))
plt.plot(mc_plt_x,np.exp(mc_a*mc_plt_x+mc_b)+mc_adj)

#for q_,_f_avg,_da_avg in zip(q_a,f_avg,da_avg):
#    if (_f_avg>0.):
#        plt.figure(3)
#        plt.scatter([len(q_)],[1./_f_avg])
#        plt.figure(4)
#        plt.scatter([_f_avg],[_da_avg])
a_starts=[]
f_avg_calc=[]
for q_,_f_avg,_da_avg,_p_info in zip(q_a,f_avg,da_avg,p_info):
    if (_f_avg>0.):
        a_starts.append(np.real(_p_info[0][1][0]))
        f_avg_calc.append(_f_avg)

# Calculate line function for thresholding initial amplitude values of partials
A_as_v_f=np.c_[np.ones(len(f_avg_calc)),f_avg_calc]
b_as_v_f=np.array(a_starts)
print A_as_v_f
print b_as_v_f
as_v_f_th=np.linalg.lstsq(A_as_v_f,b_as_v_f)[0]
plt.plot([0,max(f_avg_calc)],[as_v_f_th[0],as_v_f_th[0]+as_v_f_th[1]*max(f_avg_calc)])

for q_,_f_avg,_da_avg,_p_info in zip(q_a,f_avg,da_avg,p_info):
    if (_f_avg>0.):
        # only plot if initial amplitude above threshold
        if((as_v_f_th[0]+as_v_f_th[1]*_f_avg) < np.real(_p_info[0][1][0])): 
            plt.figure(3)
            plt.scatter(_f_avg,np.real(_p_info[0][1][0]))

plt.show()
