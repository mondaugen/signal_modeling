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

# Argument 1 is the number of frames in a group for the LP, by default 4
N_f=4
if (len(sys.argv) >= 2):
    N_f=int(sys.argv[1])

R_B=1.
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
S={}
k_S=0
F=[]
plt.figure(1)
while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,8,4,0.1,M/2))
    x_=x[h:(h+M)]
    Y[:,k]=np.fft.fft(x_*w)/w_s
    F.append(list())
    for a_ in a[-1]:
        S[k_S]=pp.LPNode(a_,[],[],k)
        F[-1].append(k_S)
        k_S+=1
        w0=np.imag(a_[1])/(2.*np.pi)*M
        w1=(np.imag(a_[1])+2.*np.imag(a_[2])*(H-1))/(2.*np.pi)*M
        n0=k
        n1=k+1.*(H-1)/H
        plt.plot([n0,n1],[w0,w1],c='k')
    k+=1
    h+=H
for k_F in xrange(len(F)-1):
    for f in F[k_F]:
        S[f].out_nodes=F[k_F+1]
    for f in F[k_F+1]:
        S[f].in_nodes=F[k_F]

def _cost_func(a,b):
    # The cost of predicting the frequency of b using the parameters of a
    return (((np.imag(a.value[1]) + np.imag(a.value[2])*H) -
            np.imag(b.value[1]))**2.) + 1.

sols=[]
sols_S=[]
sols_F=[]

for f_k in xrange(len(F)-N_f+1):
    sys.stderr.write('Solving LPs for hop %d.\n' % (f_k,))
    F_=F[f_k:f_k+N_f]
    if min([len(f_) for f_ in F_]) == 0:
        continue
    f_min=min([min(f_) for f_ in F_])
    F_=[[f__ - f_min for f__ in f_] for f_ in F_]
    S_=dict()
    for f_k_ in xrange(1,N_f-1):
        for f_ in F_[f_k_]:
            S_[f_]=deepcopy(S[f_+f_min])
            S_[f_].in_nodes=F_[f_k_-1]
            S_[f_].out_nodes=F_[f_k_+1]
            S_[f_].frame_num=f_k_
    for f_ in F_[0]:
        S_[f_]=deepcopy(S[f_+f_min])
        S_[f_].in_nodes=[]
        S_[f_].out_nodes=F_[1]
        S_[f_].frame_num=0
    for f_ in F_[N_f-1]:
        S_[f_]=deepcopy(S[f_+f_min])
        S_[f_].in_nodes=F_[N_f-2]
        S_[f_].out_nodes=[]
        S_[f_].frame_num=N_f-1
    # find the most likely connected paths using linear programming
    J_min=int(min([len(f) for f in F_])*R_B)
    if (J_min < 1):
        J_min=1
    sys.stderr.write('\tnumber of paths: %d.\n' % (J_min,))
    d=pp.g_f_2lp(S_,F_,J_min,_cost_func,{'calc_mean':0,'min_mean_dev':0})
    sols.append(solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b']))
    sols_S.append(S_)
    sols_F.append(F_)

plt.figure(2)
for k,S_,F_ in zip(xrange(len(sols)),sols_S,sols_F):
    sol=sols[k]
    N_nodes=len(S_.keys())
    for f_k in xrange(len(F_)):
        for f_ in F_[f_k]:
            for dest in S_[f_].out_nodes:
                if sol['x'][dest*N_nodes+f_] > 0.5:
                    #a_=S_[f_].value
                    #w0=np.imag(a_[1])/(2.*np.pi)*M
                    #w1=(np.imag(a_[1])+2.*np.imag(a_[2])*(H-1))/(2.*np.pi)*M
                    #n0=f_k+k
                    #n1=n0+1.*(H-1)/H
                    #plt.plot([n0,n1],[w0,w1],c='g')
                    a_0=S_[f_].value
                    a_1=S_[dest].value
                    w0_0=np.imag(a_0[1])/(2.*np.pi)*M
                    w0_1=np.imag(a_1[1])/(2.*np.pi)*M
                    n0=k+S_[f_].frame_num
                    n1=k+S_[dest].frame_num
                    plt.plot([n0,n1],[w0_0,w0_1],'g')

plt.figure(1)
plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')
plt.figure(2)
plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')
plt.show()




#
#
#for f_k in xrange(len(F)):
#    for f_ in F[f_k]:
#        for dest in S[f_].out_nodes:
#            if sol['x'][dest*N_nodes+f_] > 0.5:
#                sys.stderr.write('%f\n' % (sol['x'][dest*N_nodes+f_],))
#                a_=S[f_].value
#                w0=np.imag(a_[1])/(2.*np.pi)*M
#                w1=(np.imag(a_[1])+2.*np.imag(a_[2])*(H-1))/(2.*np.pi)*M
#                n0=f_k
#                n1=f_k+1.*(H-1)/H
#                plt.plot([n0,n1],[w0,w1],c='g')
#
