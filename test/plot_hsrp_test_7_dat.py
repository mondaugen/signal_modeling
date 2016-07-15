import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import ptpath
import ptpath_test
import cvxopt
from cvxopt import solvers
import sigmod as sm

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

#show_plots=False
show_plots=True
infilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
infilepath+='hsrp_test_7.dat'
outfilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
outfilepath+='hsrp_test_7_plot_'

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
# 1st estimated source only, best path connection
fig7=plt.figure(7)
ax7=fig7.gca()
# 2nd estimated source only, best path connection
fig8=plt.figure(8)
ax8=fig8.gca()
# 1st estimated source only, J best path connection
fig9=plt.figure(9)
ax9=fig9.gca()
# 2nd estimated source only, J best path connection
fig10=plt.figure(10)
ax10=fig10.gca()

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
ax1.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax1.set_title('Original data-points')
ax1.set_xlim(0,(h-H)/float(Fs))
ax1.legend(lgd,('Source 1','Source 2'))
fig1.savefig(outfilepath+'orig_data.eps')

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
fig2.savefig(outfilepath+'orig_spur_data.eps')

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
fig3.savefig(outfilepath+'class_data.eps')

h=0
lgd=[[],[],[]]
for i in xrange(len(datin)):
    A_pca=datin[i]['A_pca'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    m0=datin[i]['m0'][0][0]
    s0=datin[i]['s0'][0][0]
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
fig4.savefig(outfilepath+'class_pcs.eps')

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
ax5.set_title('Source 1 without swaps')
ax5.legend(lgd,('Estimated','True'))
ax5.set_xlim(0,(h-H)/float(Fs))
ax5.set_ylim(0,3.5)
fig5.savefig(outfilepath+'source_1.eps')

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
ax6.set_title('Source 2 without swaps')
ax6.legend(lgd,('Estimated','True'))
ax6.set_xlim(0,(h-H)/float(Fs))
ax6.set_ylim(0,3.5)
fig6.savefig(outfilepath+'source_2.eps')

# Arrays of start and end frequencies for each frame, each column corresponds to
# a priori source
W_ap=list()
for i in xrange(len(datin)):
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    w0_s1=[]
    w1_s1=[]
    w0_s2=[]
    w1_s2=[]
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        if (r[4]=='b'):
            w0_s1.append(float(r[2]))
            w1_s1.append(float(r[3]))
        elif (r[4]=='c'):
            w0_s2.append(float(r[2]))
            w1_s2.append(float(r[3]))
    W_ap.append([np.array([w0_s1,w1_s1]),np.array([w0_s2,w1_s2])])

# list storing the uncertainty of the classification of a set of partials
uc=list()
for i in xrange(len(W_ap)):
    # index 0 is uncertainty in keeping source classifications the way they are
    # index 1 is uncertainty in swapping them
    uc.append(np.zeros((2,),dtype='float64'))
    for c_W_ap in xrange(len(W_ap[i])):
        for c_1_W_ap in xrange(len(W_ap[i])):
            uc_=0.
            if (i>0):
                # Backwards uncertainty
                D_=np.abs(np.subtract.outer(W_ap[i][c_W_ap][0,:],W_ap[i-1][c_1_W_ap][1,:]))
                mindim=np.array(D_.shape).argmin()
                D_min=D_.min(mindim)
                uc_+=np.sum(D_min)/len(D_min)
            if (i<(len(W_ap)-1)):
                # Forwards uncertainty
                D_=np.abs(np.subtract.outer(W_ap[i][c_W_ap][1,:],W_ap[i+1][c_1_W_ap][0,:]))
                mindim=np.array(D_.shape).argmin()
                D_min=D_.min(mindim)
                uc_+=np.sum(D_min)/len(D_min)
            # if on second column, comparing with first column is swapping,
            # if on first column, comparing with second column is swapping
            swp_idx=(c_W_ap+c_1_W_ap)%2
            uc[-1][swp_idx]+=uc_

# LP to solve best J swaps
# solve
# minimize gain in uncertainty in swapping, i.e.,
# min (uc[:,1]-uc[:,0])^T*x
# s.t
# 1^T*x = J
# 0 < x < 1
J_lp=3.
uc_ary=np.array(uc)
c_lp=cvxopt.matrix(uc_ary[:,1]-uc_ary[:,0])
A_lp=cvxopt.matrix(np.array([np.ones(len(c_lp))]))
b_lp=cvxopt.matrix(np.array([[J_lp]]))
G_lp=cvxopt.matrix(np.concatenate((np.eye(len(c_lp)),-1.*np.eye(len(c_lp)))))
h_lp=cvxopt.matrix(np.concatenate((np.ones((len(c_lp),)),np.zeros((len(c_lp),)))))
sol_lp=solvers.lp(c_lp,G_lp,h_lp,A_lp,b_lp)
q_likely_lp=sol_lp['x']
#print q_likely
q_likely=[uc_.argmin() for uc_ in uc]

#Li_Ss=list()
#Li_F_s=list()
#Li_Xs=list()
#Li_sols=list()
#Li_paths=list()
#Li_Js=list()
#C=list()
#C_cxn=list()
#A_vit=list()
#def _node_cxn_cost(a,b,dt):
#    return ((a[0]+a[1]*dt)-b[0])**2.
#def _lp_cost_fun(a,b):
#    return _node_cxn_cost(a.value,b.value,(b.frame_num-a.frame_num)*H/float(Fs))+1.
#h=0
#sols=[]
#lgd=[[],[]]
#for i in xrange(len(datin)-1):
#    X_plt_est=datin[i]['X_plt'][0][0]
#    clr_e=datin[i]['clr_'][0][0]
#    X_s1_0=[]
#    X_s2_0=[]
#    X_s1_1=[]
#    X_s2_1=[]
#    for r,s in zip(X_plt_est,clr_e):
#        if (s=='b'):
#            X_s1_0.append(np.array([r[0],r[1]]))
#        if (s=='c'):
#            X_s2_0.append(np.array([r[0],r[1]]))
#    X_plt_est=datin[i+1]['X_plt'][0][0]
#    clr_e=datin[i+1]['clr_'][0][0]
#    for r,s in zip(X_plt_est,clr_e):
#        if (s=='b'):
#            X_s1_1.append(np.array([r[0],r[1]]))
#        if (s=='c'):
#            X_s2_1.append(np.array([r[0],r[1]]))
#
#    # LPs, find best connections
#    # 0=1-1 1=2-2 2=1-2 3=2-1
#    Ss=[dict() for _ in xrange(4)]
#    F_s=[[] for _ in xrange(4)]
#    Xs=[(X_s1_0,X_s1_1),
#        (X_s2_0,X_s2_1),
#        (X_s1_0,X_s2_1),
#        (X_s2_0,X_s1_1)]
#    sols=[]
#    paths=[]
#    Js=[]
#    for S,F,X in zip(Ss,F_s,Xs):
#        n_node=0
#        n_frame=0
#        # Build frames and graph
#        for x_ in X:
#            F.append([])
#            for k_a_ in xrange(len(x_)):
#                x__=x_[k_a_]
#                F[-1].append(n_node)
#                S[n_node]=ptpath.LPNode(value=x__,frame_num=n_frame)
#                n_node+=1
#            n_frame+=1
#        for f in F[0]:
#            S[f].out_nodes=F[1]
#        for f in F[1]:
#            S[f].in_nodes=F[0]
#        J=min(len(X_) for X_ in X)
#        Js.append(J)
#        # Build linear program
#        d=ptpath.g_f_2lp(S,F,J,_lp_cost_fun,{'calc_mean':0,'min_mean_dev':0})
#        # Solve LP
#        sols.append(solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b']))
#        # Extract paths
#        paths.append(ptpath_test.lp_sol_extract_paths(sols[-1]['x'],S,F))
#    # add to list
#    Li_Ss.append(Ss)
#    Li_F_s.append(F_s)
#    Li_Xs.append(Xs)
#    Li_sols.append(sols)
#    Li_paths.append(paths)
#    Li_Js.append(Js)
#    # Add costs to cost lattice
#    C.append(np.array([sols[0]['primal objective']/Js[0]
#                +sols[1]['primal objective']/Js[1],
#                 sols[2]['primal objective']/Js[2]
#                +sols[3]['primal objective']/Js[3]]))
#    # Add connections to connections lattice
#    C_cxn.append([((0,1),(0,1)),((0,1),(1,0))])
#    # Add transition cost matrix to A_vit
#    A_vit.append(np.array([
#        [sols[0]['primal objective'],sols[2]['primal objective']],
#        [sols[3]['primal objective'],sols[0]['primal objective']]
#        ]))

#C=list()
#for i in xrange(len(datin)-1):
#    X_plt_est=datin[i]['X_plt'][0][0]
#    clr_e=datin[i]['clr_'][0][0]
#    X_s1_0=[]
#    X_s2_0=[]
#    X_s1_1=[]
#    X_s2_1=[]
#    for r,s in zip(X_plt_est,clr_e):
#        if (s=='b'):
#            X_s1_0.append(np.array([r[0],r[1]]))
#        if (s=='c'):
#            X_s2_0.append(np.array([r[0],r[1]]))
#    X_plt_est=datin[i+1]['X_plt'][0][0]
#    clr_e=datin[i+1]['clr_'][0][0]
#    for r,s in zip(X_plt_est,clr_e):
#        if (s=='b'):
#            X_s1_1.append(np.array([r[0],r[1]]))
#        if (s=='c'):
#            X_s2_1.append(np.array([r[0],r[1]]))
#    #costs=np.ndarray((2,),dtype='float64')
#    #Xs=[[(np.array(X_s1_0),np.array(X_s1_1)),
#    #     (np.array(X_s2_0),np.array(X_s2_1))],
#    #    [(np.array(X_s1_0),np.array(X_s2_1)),
#    #     (np.array(X_s2_0),np.array(X_s1_1))]]
#    #for xs,costs_i in zip(Xs,xrange(len(costs))):
#    #    cost=0.
#    #    ncxns=0
#    #    for xt0,xt1 in xs:
#    #        # Carry out greedy search for shortest paths
#    #        D=np.subtract.outer(xt0[:,0]+xt0[:,1]*H,xt1[:,0])**2.
#    #        r=list(xrange(D.shape[0]))
#    #        c=list(xrange(D.shape[1]))
#    #        while (True):
#    #            min_rc=(-1,-1)
#    #            min_cost=float('inf')
#    #            for r_ in r:
#    #                for c_ in c:
#    #                    if D[r_,c_] < min_cost:
#    #                        min_cost=D[r_,c_]
#    #                        min_rc=(r_,c_)
#    #            if (min_rc[0]==-1):
#    #                break
#    #            cost+=min_cost
#    #            ncxns+=1
#    #            r.remove(min_rc[0])
#    #            c.remove(min_rc[1])
#    #    costs[costs_i]=cost/float(ncxns)
#    costs=np.ndarray((2,2),dtype='float64')
#    Xs=[[(np.array(X_s1_0),np.array(X_s1_1)),(np.array(X_s2_0),np.array(X_s1_1))],
#        [(np.array(X_s1_0),np.array(X_s2_1)),(np.array(X_s2_0),np.array(X_s2_1))]]
#    for r_xs in xrange(2):
#        for c_xs in xrange(2):
#            xt0,xt1=Xs[r_xs][c_xs]
#            cost=0.
#            ncxns=0
#            # Carry out greedy search for shortest paths
#            D=np.subtract.outer(xt0[:,0]+xt0[:,1]*H,xt1[:,0])**2.
#            r=list(xrange(D.shape[0]))
#            c=list(xrange(D.shape[1]))
#            while (True):
#                min_rc=(-1,-1)
#                min_cost=float('inf')
#                for r_ in r:
#                    for c_ in c:
#                        if D[r_,c_] < min_cost:
#                            min_cost=D[r_,c_]
#                            min_rc=(r_,c_)
#                if (min_rc[0]==-1):
#                    break
#                cost+=min_cost
#                ncxns+=1
#                r.remove(min_rc[0])
#                c.remove(min_rc[1])
#            costs[r_xs,c_xs]=cost/float(ncxns)
#    C.append(costs)

# Now carry out viterbi algorithm to find the 2 most likely paths
# (This could be done with the classic viterbi algorithm but it's easier to set
# up this way)
#q_likely=ptpath_test.shortest_paths_viterbi(C,C_cxn)
#q_likely,p_likely=sm.best_path_viterbi(C)
#q_likely=np.concatenate(([0],q_likely))
#q_likely=[c_.argmin() for c_ in C]
#q_likely=np.concatenate(([0],q_likely))
#q_likely=np.zeros((32,))
#q_likely[:16]=np.array([0,0,0,0,0,1,1,0,0,0,1,1,0,0,1,1])
#q_likely[:15]=np.array([0,0,0,0,1,1,0,0,0,1,1,0,0,1,1])
#q_likely=[uc_.argmin() for uc_ in uc]

#h=0
#lgd=[[],[]]
#_plt_clr='b'
#_src=0
#for i,q_,lip_,lis_ in zip(xrange(len(datin)-1),
#                          q_likely[:-1],
#                          Li_paths[:-1],
#                          Li_Ss[:-1]):
#    # swap if indicated by best path
#    if (q_==1):
#        # Source 1 'b', source 2 'c'
#        if (_plt_clr=='b'):
#            _plt_clr='c'
#            _src=2
#        else:
#            _plt_clr='b'
#            _src=0
#    w0_e=[]
#    w1_e=[]
#    for pth in lip_[_src]:
#        if (len(pth) < 2):
#            continue
#        w0_e.append(lis_[_src][pth[0]].value[0])
#        w1_e.append(lis_[_src][pth[1]].value[0])
#    w0_e=np.array(w0_e)
#    w1_e=np.array(w1_e)
#    N_w0_e=len(w0_e)
#    h_e=h+H*np.c_[np.zeros(N_w0_e),np.ones(N_w0_e)]
#    h_e/=float(Fs)
#    for r in np.c_[h_e,w0_e,w1_e]:
#        ls_=ls_s1
#        lgd_idx=0
#        lgd[lgd_idx]=ax7.plot(r[0:2],r[2:4],ls_)[0]
#    X_plt_orig=datin[i]['X_plt_orig'][0][0]
#    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
#    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
#    N_w0_o=len(w0_o)
#    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
#    h_o/=Fs
#    lgd[1]=ax7.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
#    h+=H
#ax7.set_xlabel('Time (seconds)')
#ax7.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#ax7.set_title('Source 1 after connecting paths')
#ax7.legend(lgd,('Estimated','True'))
#ax7.set_xlim(0,(h-H)/float(Fs))
#ax7.set_ylim(0,3.5)

#h=0
#lgd=[[]]
#_clr='b'
#for i,q_ in zip(xrange(len(datin)),q_likely):
#    X_plt_est=datin[i]['X_plt'][0][0]
#    clr_e=datin[i]['clr_'][0][0]
#    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
#    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
#    N_w0_e=len(w0_e)
#    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
#    h_e/=float(Fs)
#    if (q_==1):
#        if (_clr == 'b'):
#            _clr='c'
#        else:
#            _clr='b'
#    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
#        if (r[4]==_clr):
#            ls_=ls_s1
#            lgd_idx=0
#            lgd[lgd_idx]=ax7.plot(r[0:2],r[2:4],ls_)[0]
#    X_plt_orig=datin[i]['X_plt_orig'][0][0]
#    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
#    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
#    N_w0_o=len(w0_o)
#    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
#    h_o/=Fs
##    lgd[1]=ax7.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
#    h+=H
#ax7.set_xlabel('Time (seconds)')
#ax7.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#ax7.set_title('Source 1 After Conditionally Swapping Sources')
##ax7.legend(lgd,('Estimated','True'))
#ax7.legend(lgd,('Estimated',))
#ax7.set_xlim(0,(h-H)/float(Fs))
#ax7.set_ylim(0,3.5)

h=0
lgd=[[],[]]
for i,q_ in zip(xrange(len(datin)),q_likely):
    q_=int(round(q_))
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    if (q_==1):
        _clr='c'
    else:
        _clr='b'
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        if (r[4]==_clr):
            ls_=ls_s1
            lgd_idx=0
            lgd[lgd_idx]=ax7.plot(r[0:2],r[2:4],ls_)[0]
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    h_o/=Fs
    lgd[1]=ax7.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
    h+=H
ax7.set_xlabel('Time (seconds)')
ax7.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax7.set_title('Source 1 after all gainful swaps')
ax7.legend(lgd,('Estimated','True'))
#ax7.legend(lgd,('Estimated',))
ax7.set_xlim(0,(h-H)/float(Fs))
ax7.set_ylim(0,3.5)
fig7.savefig(outfilepath+'source_1_all_swapped.eps')

h=0
lgd=[[],[]]
for i,q_ in zip(xrange(len(datin)),q_likely):
    q_=int(round(q_))
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    if (q_==1):
        _clr='b'
    else:
        _clr='c'
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        if (r[4]==_clr):
            ls_=ls_s2
            lgd_idx=0
            lgd[lgd_idx]=ax8.plot(r[0:2],r[2:4],ls_)[0]
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    h_o/=Fs
    lgd[1]=ax8.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
    h+=H
ax8.set_xlabel('Time (seconds)')
ax8.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax8.set_title('Source 2 after all gainful swaps')
ax8.legend(lgd,('Estimated','True'))
#ax8.legend(lgd,('Estimated',))
ax8.set_xlim(0,(h-H)/float(Fs))
ax8.set_ylim(0,3.5)
fig8.savefig(outfilepath+'source_2_all_swapped.eps')

h=0
lgd=[[],[]]
for i,q_ in zip(xrange(len(datin)),q_likely_lp):
    q_=int(round(q_))
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    if (q_==1):
        _clr='c'
    else:
        _clr='b'
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        if (r[4]==_clr):
            ls_=ls_s1
            lgd_idx=0
            lgd[lgd_idx]=ax9.plot(r[0:2],r[2:4],ls_)[0]
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    h_o/=Fs
    lgd[1]=ax9.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
    h+=H
ax9.set_xlabel('Time (seconds)')
ax9.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax9.set_title('Source 1 after %d most gainful swaps' % (int(J_lp),))
ax9.legend(lgd,('Estimated','True'))
#ax9.legend(lgd,('Estimated',))
ax9.set_xlim(0,(h-H)/float(Fs))
ax9.set_ylim(0,3.5)
fig9.savefig(outfilepath+'source_1_some_swapped.eps')

h=0
lgd=[[],[]]
for i,q_ in zip(xrange(len(datin)),q_likely_lp):
    q_=int(round(q_))
    X_plt_est=datin[i]['X_plt'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
    w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
    N_w0_e=len(w0_e)
    h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
    h_e/=float(Fs)
    if (q_==1):
        _clr='b'
    else:
        _clr='c'
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        if (r[4]==_clr):
            ls_=ls_s2
            lgd_idx=0
            lgd[lgd_idx]=ax10.plot(r[0:2],r[2:4],ls_)[0]
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    h_o/=Fs
    lgd[1]=ax10.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
    h+=H
ax10.set_xlabel('Time (seconds)')
ax10.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax10.set_title('Source 2 after %d most gainful swaps' % (int(J_lp,)))
ax10.legend(lgd,('Estimated','True'))
ax10.set_xlim(0,(h-H)/float(Fs))
ax10.set_ylim(0,3.5)
fig10.savefig(outfilepath+'source_2_some_swapped.eps')

if (show_plots):
    plt.show()
