import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import ptpath
import ptpath_test
from cvxopt import solvers

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

show_plots=False
infilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
infilepath+='hsrp_test_7.dat'

# Original data
fig1=plt.figure(1)
ax1=fig1.gca()
# Original + spurious
fig2=plt.figure(2)
ax2=fig2.gca()
# Original + spurious, estimated membership
fig3=plt.figure(3)
ax3=fig3.gca()
# 1st PC
fig4=plt.figure(4)
ax4=fig4.gca()
# 1st estimated source only
fig5=plt.figure(5)
ax5=fig5.gca()
# 2nd estimated source only
fig6=plt.figure(6)
ax6=fig6.gca()

# Hop size
H=256
Fs=16000
K=20

# Line styles
# source 1
ls_s1='k:'
# source 2
ls_s2='k--'
# background
ls_bg='lightgrey'
# Point styles
ps_s1_s='o'
ps_s1_c='k'
ps_s2_s='o'
ps_s2_c='grey'
ps_bg_c='lightgrey'
ps_bg_s='.'

datin=sio.loadmat(infilepath)['datout'][0]
h=0
lgd=[[],[]]
for i in xrange(len(datin)):
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    h_o/=Fs
    lgd[0]=ax1.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_s1)[0]
    lgd[1]=ax1.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_s2)[0]
    h+=H
ax1.set_xlabel('Time (seconds)')
# Why undefined sequence?
ax1.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax1.set_title('Original data-points')
ax1.set_xlim(0,(h-H)/float(Fs))
ax1.legend(lgd,('Source 1','Source 2'))

h=0
for i in xrange(len(datin)):
    X_plt=datin[i]['X_plt'][0][0]
    w0=X_plt[:,0]+(-0.5)*X_plt[:,1]*H
    w1=X_plt[:,0]+(0.5)*X_plt[:,1]*H
    N_w0=len(w0)
    h_=h+H*np.c_[-0.5*np.ones(N_w0),0.5*np.ones(N_w0)]
    h_/=float(Fs)
    ax2.plot(h_.T,np.c_[w0,w1].T,'k')
    h+=H
ax2.set_xlim(0,(h-H)/float(Fs))
ax2.set_ylim(0,3.5)
ax2.set_xlabel('Time (seconds)')
ax2.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax2.set_title('Original and spurious data-points')

h=0
lgd=[[],[],[]]
for i in xrange(len(datin)):
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        # Source 1 'b'
        if (r[4]=='b'):
            ls_=ls_s1
            lgd_idx=0
        # Source 2 'c'
        elif (r[4]=='c'):
            ls_=ls_s2
            lgd_idx=1
        else:
            ls_=ls_bg
            lgd_idx=2
        lgd[lgd_idx]=ax3.plot(r[0:2],r[2:4],ls_)[0]
    h+=H
ax3.set_xlabel('Time (seconds)')
ax3.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax3.set_title('Classified data-points')
ax3.legend(lgd,('Source 1','Source 2','Spurious'))
ax3.set_xlim(0,(h-H)/float(Fs))
ax3.set_ylim(0,3.5)

h=0
lgd=[[],[],[]]
for i in xrange(len(datin)):
    A_pca=datin[i]['A_pca'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    a0=A_pca[:,0]
    N_a0=len(a0)
    h_e=h+H*np.c_[np.ones(N_a0)]
    h_e/=float(Fs)
    n_s1=0
    n_s2=0
    for r in np.c_[h_e,a0,clr_e]:
        # Source 1 'b'
        if (r[2]=='b'):
            ls_s=ps_s1_s
            ls_c=ps_s1_c
            lw_=0
            lgd_idx=0
            n_s1+=1
        # Source 2 'c'
        elif (r[2]=='c'):
            ls_s=ps_s2_s
            ls_c=ps_s2_c
            lw_=0
            lgd_idx=1
            n_s2+=1
        else:
            ls_s=ps_bg_s
            ls_c=ps_bg_c
            lw_=0
            lgd_idx=2
        lgd[lgd_idx]=ax4.scatter(r[0],r[1],marker=ls_s,c=ls_c,lw=lw_)
    print n_s1,n_s2
    h+=H
ax4.set_xlabel('Time (seconds)')
ax4.set_ylabel('1st PC')
ax4.set_title('Principal components and their classification')
ax4.legend(lgd,('Source 1','Source 2','Spurious'))
ax4.set_xlim(0.2,0.3)
ax4.set_ylim(-0.0012,0.0012)

h=0
lgd=[[],[]]
for i in xrange(len(datin)):
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        # Source 1 'b'
        if (r[4]=='b'):
            ls_=ls_s1
            lgd_idx=0
            lgd[lgd_idx]=ax5.plot(r[0:2],r[2:4],ls_)[0]
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    h_o/=Fs
    lgd[1]=ax5.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
    h+=H
ax5.set_xlabel('Time (seconds)')
ax5.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax5.set_title('Source 1 before connecting paths')
ax5.legend(lgd,('Estimated','True'))
ax5.set_xlim(0,(h-H)/float(Fs))
ax5.set_ylim(0,3.5)

h=0
lgd=[[],[]]
for i in xrange(len(datin)):
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        # Source 2 'c'
        if (r[4]=='c'):
            ls_=ls_s2
            lgd_idx=0
            lgd[lgd_idx]=ax6.plot(r[0:2],r[2:4],ls_)[0]
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    h_o/=Fs
    lgd[1]=ax6.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
    h+=H
ax6.set_xlabel('Time (seconds)')
ax6.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax6.set_title('Source 2 before connecting paths')
ax6.legend(lgd,('Estimated','True'))
ax6.set_xlim(0,(h-H)/float(Fs))
ax6.set_ylim(0,3.5)

Li_Ss=list()
Li_F_s=list()
Li_Xs=list()
Li_sols=list()
Li_paths=list()
C=list()
C_cxn=list()
def _node_cxn_cost(a,b,dt):
    return ((a[0]+a[1]*dt)-b[0])**2.
def _lp_cost_fun(a,b):
    return _node_cxn_cost(a.value,b.value,(b.frame_num-a.frame_num)*H/float(Fs))+1.
h=0
sols=[]
lgd=[[],[]]
for i in xrange(len(datin)-1):
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    X_s1_0=[]
    X_s2_0=[]
    X_s1_1=[]
    X_s2_1=[]
    for r,s in zip(X_plt_est,clr_e):
        if (s=='b'):
            X_s1_0.append(np.array([r[0],r[1]]))
        if (s=='c'):
            X_s2_0.append(np.array([r[0],r[1]]))
    X_plt_est=datin[i+1]['X_plt'][0][0]
    clr_e=datin[i+1]['clr_'][0][0]
    for r,s in zip(X_plt_est,clr_e):
        if (s=='b'):
            X_s1_1.append(np.array([r[0],r[1]]))
        if (s=='c'):
            X_s2_1.append(np.array([r[0],r[1]]))

    # LPs, find best connections
    # 0=1-1 1=2-2 2=1-2 3=2-1
    Ss=[dict() for _ in xrange(4)]
    F_s=[[] for _ in xrange(4)]
    Xs=[(X_s1_0,X_s1_1),
        (X_s2_0,X_s2_1),
        (X_s1_0,X_s2_1),
        (X_s2_0,X_s1_1)]
    sols=[]
    paths=[]
    for S,F,X in zip(Ss,F_s,Xs):
        n_node=0
        n_frame=0
        # Build frames and graph
        for x_ in X:
            F.append([])
            for k_a_ in xrange(len(x_)):
                x__=x_[k_a_]
                F[-1].append(n_node)
                S[n_node]=ptpath.LPNode(value=x__,frame_num=n_frame)
                n_node+=1
            n_frame+=1
        for f in F[0]:
            S[f].out_nodes=F[1]
        for f in F[1]:
            S[f].in_nodes=F[0]
        J=min(len(X_) for X_ in X)
        # Build linear program
        d=ptpath.g_f_2lp(S,F,J,_lp_cost_fun,{'calc_mean':0,'min_mean_dev':0})
        # Solve LP
        sols.append(solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b']))
        # Extract paths
        paths.append(ptpath_test.lp_sol_extract_paths(sols[-1]['x'],S,F))
    # add to list
    Li_Ss.append(Ss)
    Li_F_s.append(F_s)
    Li_Xs.append(Xs)
    Li_sols.append(sols)
    Li_paths.append(paths)
    # Add costs to cost lattice
    C.append([sols[0]['primal objective']+sols[1]['primal objective'],
              sols[2]['primal objective']+sols[3]['primal objective']])
    # Add connections to connections lattice
    C_cxn.append([((0,1),(0,1)),((0,1),(1,0))])

# Now carry out viterbi algorithm to find the 2 most likely paths
# (This could be done with the classic viterbi algorithm but it's easier to set
# up this way)
q_likely=ptpath_test.shortest_paths_viterbi(C,C_cxn)

if (show_plots):
    plt.show()
