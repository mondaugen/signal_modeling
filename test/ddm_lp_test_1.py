import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import sys
import ptpath as pp
from cvxopt import solvers

H=256
#N=4096
x=np.fromfile(sys.stdin)
N=len(x)
n=np.arange(N)
#plt.plot(n,np.real(x))
a=[]
h=0
M=1024
Y=np.zeros((M,(N-M)/H+1)).astype('complex_')
m=np.arange(M)
w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
k=0
w_s=np.sum(w)
S={}
k_S=0
F=[]
while ((h+M) <= N):
    sys.stderr.write('%d / %d\n' % (h,N))
    a.append(sm.ddm_p2_1_3_b(x[h:(h+M)],w,dw,6,4,0.01,M/2))
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

# find the most likely connected paths using linear programming
J_min=min([len(f) for f in F])
d=pp.g_f_2lp(S,F,J_min,_cost_func,{'calc_mean':0,'min_mean_dev':0})
sol=solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b'])
N_nodes=len(S.keys())

for f_k in xrange(len(F)):
    for f_ in F[f_k]:
        for dest in S[f_].out_nodes:
            if sol['x'][dest*N_nodes+f_] > 0.5:
                sys.stderr.write('%f\n' % (sol['x'][dest*N_nodes+f_],))
                a_=S[f_].value
                w0=np.imag(a_[1])/(2.*np.pi)*M
                w1=(np.imag(a_[1])+2.*np.imag(a_[2])*(H-1))/(2.*np.pi)*M
                n0=f_k
                n1=f_k+1.*(H-1)/H
                plt.plot([n0,n1],[w0,w1],c='g')

plt.imshow(20.*np.log10(np.abs(Y)),origin='lower',aspect='auto',interpolation='nearest')
plt.show()
