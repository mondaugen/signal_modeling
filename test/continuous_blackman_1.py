import cvxopt
import numpy as np
import matplotlib.pyplot as plt
import os

show_plots=False

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

foutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
P=cvxopt.matrix(np.eye(4))
# Minimum 4-term blackman-harris
a_=np.array([.35875,.48829,.14128,.01168])
#a_=np.array([.40217,.49703,.09392,.00183])
q=-1.*cvxopt.matrix(a_)
A=cvxopt.matrix([[1,-1,1,-1],[1,1,1,1]],tc='d').T
b=cvxopt.matrix([0,1],tc='d')
# Continuous blackman-harris based on min 4-term
a=cvxopt.solvers.qp(P,q,A=A,b=b)['x']
N=32
n=np.arange(N)
a=np.array(a)
print a.flatten()
w=np.sum(np.cos(np.outer(2.*np.pi*n/N,[0,1,2,3]))*np.power(-1,[0,1,2,3])*a.flatten(),
        axis=1)
w2=np.sum(np.cos(np.outer(2.*np.pi*n/N,[0,1,2,3]))*np.power(-1,[0,1,2,3])*a_.flatten(),
       axis=1)
print w
print w[0],w[N-1]
ov=64
W=np.fft.fft(w,ov*N)/sum(w)
W_orig=W
W2=np.fft.fft(w2,ov*N)/sum(w2)
n-=N/2
n__=np.linspace(0,N,N*ov)
w_=np.sum(np.cos(np.outer(2.*np.pi*n__/N,[0,1,2,3]))*np.power(-1,[0,1,2,3])*a.flatten(),
        axis=1)
w2_=np.sum(np.cos(np.outer(2.*np.pi*n__/N,[0,1,2,3]))*np.power(-1,[0,1,2,3])*a_.flatten(),
       axis=1)
n__-=N/2
n_=np.arange(ov*N)/float(ov)
n_-=N/2
_N=len(W)
W=np.r_[W[_N/2:],W[:_N/2]]
W2=np.r_[W2[_N/2:],W2[:_N/2]]
plt.figure(1)
plt.plot(n__,w_,c='k')
plt.title('$\mathcal{C}^{1}$ 4-Term Blackman-Harris: Time Domain')
plt.xlabel('Sample number')
plt.ylabel('Sample value')
ax1=plt.gca()
ax1.set_ylim(0,1.1)
ax1.set_xlim(-N/2,N/2)
plt.savefig(foutpath+'c1_blackman_td.eps')
plt.figure(2)
plt.plot(n_/float(N)*2.*np.pi,20*np.log10(np.abs(W)),c='k')
plt.title('$\mathcal{C}^{1}$ 4-Term Blackman-Harris: Frequency Domain')
plt.xlabel('Frequency (Radians/Second)')
plt.ylabel('Power (dB)')
ax2=plt.gca()
ax2.set_ylim(-140,5)
ax2.set_xlim(-np.pi,np.pi)
plt.savefig(foutpath+'c1_blackman_fd.eps')
plt.figure(3)
plt.plot(n__,w2_,c='k')
plt.title('Minimum 4-Term Blackman-Harris: Time Domain')
plt.xlabel('Sample number')
plt.ylabel('Sample value')
ax3=plt.gca()
ax3.set_ylim(0,1.1)
ax3.set_xlim(-N/2,N/2)
plt.savefig(foutpath+'min4_blackman_td.eps')
plt.figure(4)
plt.plot(n_/float(N)*2.*np.pi,20*np.log10(np.abs(W2)),c='k')
plt.title('Minimum 4-Term Blackman-Harris: Frequency Domain')
plt.xlabel('Frequency (Radians/Second)')
plt.ylabel('Power (dB)')
ax4=plt.gca()
ax4.set_ylim(-140,5)
ax4.set_xlim(-np.pi,np.pi)
plt.savefig(foutpath+'min4_blackman_fd.eps')
plt.figure(5)
plt.plot(n_,20*np.log10(np.abs(W)),c='k')
plt.plot([min(n_),max(n_)],[-6,-6])
plt.title('Continuous blackman for examining')
print 'Height of highest side lobe for C1 window: %f' % (20*np.log10(
    np.abs(W_orig[ov*6.5])),)
print '6-dB bandwidth is 2.66 (from visual inspection)'
print 'Continuous blackman endpoints %f,%f'%(w_[0],w_[-1])
print 'Minimum blackman endpoints %f,%f' %(w2_[0],w2_[-1])
plt.figure(6)
n_cu=40
plt.plot(n__[:n_cu],w_[:n_cu],
        label='$\mathcal{C}^{1}$ 4-Term Blackman-Harris',c='grey')
plt.plot(n__[:n_cu],w2_[:n_cu],label='Minimum 4-Term Blackman-Harris',c='k')
plt.legend(loc='best')
plt.ylim(0,max(w_[:n_cu])*1.25)
plt.xlim(n__[0],n__[n_cu-1])
plt.title('Comparison of endpoints of window in time-domain')
plt.xlabel('Sample number')
plt.ylabel('Sample value')
plt.savefig(foutpath+'c1_vs_min_blackman_closeup.eps')

if (show_plots):
    plt.show()
