import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import ptpath
import ptpath_test
import cvxopt
from cvxopt import solvers
import sigmod as sm
import sys


plt.rc('text',usetex=True)
plt.rc('font',family='serif')
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

show_plots=False
#show_plots=True
infilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'

if (len(sys.argv) != 2):
    raise Exception('Specify file.')
if not sys.argv[1].endswith('.dat'):
    raise Exception('File must end with .dat')

infilepath+=sys.argv[1]
outfilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
outfilepath+=sys.argv[1][:sys.argv[1].rfind('.dat')]+'_plot_'
fname_idx0=0
fname_idx1=0


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
if (len(sys.argv) == 3):
    ls_s1='k-'
else:
    ls_s1='k:'
# source 2
if (len(sys.argv) == 3):
    ls_s2='k-'
else:
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
ps_m1_s='x'
ps_m1_c='k'
ps_m2_s='x'
ps_m2_c='grey'
ps_v1_s='_'
ps_v1_c='k'
ps_v2_s='_'
ps_v2_c='grey'

lgd_ax1=[[],[]]
lgd_ax3=[[],[],[]]
lgd_ax4=[[],[],[],[],[],[],[]]
lgd_ax5=[[],[]]
lgd_ax6=[[],[]]
lgd_ax7=[[],[]]
lgd_ax8=[[],[]]
lgd_ax9=[[],[]]
lgd_ax10=[[],[]]

for fname_idx_ in xrange(fname_idx0,fname_idx1+1):

    if (len(sys.argv) == 3):
        datin=sio.loadmat(infilepath % (fname_idx_,))['datout'][0]
    else:
        datin=sio.loadmat(infilepath)['datout'][0]

    h=0
    for i in xrange(len(datin)):
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax1[0]=ax1.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_s1)[0]
        lgd_ax1[1]=ax1.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_s2)[0]
        h+=H

    
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
    
    h=0
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
                lgd_ax3_idx=0
            # Source 2 'c'
            elif (r[4]=='c'):
                ls_=ls_s2
                lgd_ax3_idx=1
            else:
                ls_=ls_bg
                lgd_ax3_idx=2
            lgd_ax3[lgd_ax3_idx]=ax3.plot(r[0:2],r[2:4],ls_)[0]
        h+=H
    
    h=0
    for i in xrange(len(datin)):
        A_pca=datin[i]['A_pca'][0][0]
        clr_e=datin[i]['clr_'][0][0]
        m0=datin[i]['m0'][0][0]
        s0=datin[i]['s0'][0][0][0].T
        a0=A_pca[:,0]
        N_a0=len(a0)
        h_e=h+H*np.c_[np.ones(N_a0)]
        h_e/=float(Fs)
        n_s1=0
        n_s2=0
        # Plot initial mean estimate
        ls_s=ps_m1_s
        ls_c=ps_m1_c
        lgd_ax4[3]=ax4.scatter(h_e[0],m0[0],marker=ls_s,c=ls_c,s=80)
        ls_s=ps_m2_s
        ls_c=ps_m2_c
        lgd_ax4[4]=ax4.scatter(h_e[0],m0[1],marker=ls_s,c=ls_c,s=80)
        # Plot 1 standard deviation away from initial mean estimate using initial
        # variance estimate
        ls_s=ps_v1_s
        ls_c=ps_v1_c
        lgd_ax4[5]=ax4.scatter(h_e[0],m0[0]+np.sqrt(s0[0]),marker=ls_s,c=ls_c,s=80)
        lgd_ax4[5]=ax4.scatter(h_e[0],m0[0]-np.sqrt(s0[0]),marker=ls_s,c=ls_c,s=80)
        ls_s=ps_v2_s
        ls_c=ps_v2_c
        lgd_ax4[6]=ax4.scatter(h_e[0],m0[1]+np.sqrt(s0[1]),marker=ls_s,c=ls_c,s=80)
        lgd_ax4[6]=ax4.scatter(h_e[0],m0[1]-np.sqrt(s0[1]),marker=ls_s,c=ls_c,s=80)
        for r in np.c_[h_e,a0,clr_e]:
            # Source 1 'b'
            if (r[2]=='b'):
                ls_s=ps_s1_s
                ls_c=ps_s1_c
                lw_=0
                lgd_ax4_idx=0
                n_s1+=1
            # Source 2 'c'
            elif (r[2]=='c'):
                ls_s=ps_s2_s
                ls_c=ps_s2_c
                lw_=0
                lgd_ax4_idx=1
                n_s2+=1
            else:
                ls_s=ps_bg_s
                ls_c=ps_bg_c
                lw_=0
                lgd_ax4_idx=2
            lgd_ax4[lgd_ax4_idx]=ax4.scatter(r[0],r[1],marker=ls_s,c=ls_c,lw=lw_)
        print n_s1,n_s2
        h+=H
    
    h=0
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
                lgd_ax5_idx=0
                lgd_ax5[lgd_ax5_idx]=ax5.plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax5[1]=ax5.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
        h+=H
    
    h=0
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
                lgd_ax6_idx=0
                lgd_ax6[lgd_ax6_idx]=ax6.plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax6[1]=ax6.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
        h+=H
    
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
    J_lp=4.
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
    
    h=0
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
                lgd_ax7_idx=0
                lgd_ax7[lgd_ax7_idx]=ax7.plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax7[1]=ax7.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
        h+=H
    
    h=0
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
                lgd_ax8_idx=0
                lgd_ax8[lgd_ax8_idx]=ax8.plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax8[1]=ax8.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
        h+=H
    
    h=0
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
                lgd_ax9_idx=0
                lgd_ax9[lgd_ax9_idx]=ax9.plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax9[1]=ax9.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
        h+=H
    
    h=0
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
                lgd_ax10_idx=0
                lgd_ax10[lgd_ax10_idx]=ax10.plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax10[1]=ax10.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
        h+=H

ax6.set_xlabel('Time (seconds)')
ax6.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax6.set_title('Source 2 without swaps')
ax6.set_xlim(0,(h-H)/float(Fs))
ax6.set_ylim(0,3.5)

ax5.set_xlabel('Time (seconds)')
ax5.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax5.set_title('Source 1 without swaps')
ax5.set_xlim(0,(h-H)/float(Fs))
ax5.set_ylim(0,3.5)

ax4.set_xlabel('Time (seconds)')
ax4.set_ylabel('1st PC')
ax4.set_title('Principal components and their classification')
ax4.set_xlim(0.2,0.3)
ax4.set_ylim(-0.00035,0.)

ax3.set_xlabel('Time (seconds)')
ax3.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax3.set_title('Classified data-points')
ax3.set_xlim(0,(h-H)/float(Fs))
ax3.set_ylim(0,3.5)

ax2.set_xlim(0,(h-H)/float(Fs))
ax2.set_ylim(0,3.5)
ax2.set_xlabel('Time (seconds)')
ax2.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax2.set_title('Original and spurious data-points')

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax1.set_title('Original data-points')
ax1.set_xlim(0,(h-H)/float(Fs))

ax7.set_xlabel('Time (seconds)')
ax7.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax7.set_title('Source 1 after all gainful swaps')
ax7.set_xlim(0,(h-H)/float(Fs))
ax7.set_ylim(0,3.5)

ax8.set_xlabel('Time (seconds)')
ax8.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax8.set_title('Source 2 after all gainful swaps')
ax8.set_xlim(0,(h-H)/float(Fs))
ax8.set_ylim(0,3.5)

ax9.set_xlabel('Time (seconds)')
ax9.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax9.set_title('Source 1 after %d most gainful swaps' % (int(J_lp),))
ax9.set_xlim(0,(h-H)/float(Fs))
ax9.set_ylim(0,3.5)

ax10.set_xlabel('Time (seconds)')
ax10.set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
ax10.set_title('Source 2 after %d most gainful swaps' % (int(J_lp,)))
ax10.set_xlim(0,(h-H)/float(Fs))
ax10.set_ylim(0,3.5)

ax1.legend(lgd_ax1,('Source 1','Source 2'))
ax3.legend(lgd_ax3,('Source 1','Source 2','Spurious'))
ax4.legend(lgd_ax4,('Source 1','Source 2','Spurious','$\\mu_{0}^{0}$',
                    '$\\mu_{1}^{0}$','$\\pm \\sigma_{0}^{0}$',
                    '$\\pm \\sigma_{1}^{0}$'))
ax5.legend(lgd_ax5,('Estimated','True'))
ax6.legend(lgd_ax6,('Estimated','True'))
ax7.legend(lgd_ax7,('Estimated','True'))
ax8.legend(lgd_ax8,('Estimated','True'))
ax9.legend(lgd_ax9,('Estimated','True'))
ax10.legend(lgd_ax10,('Estimated','True'))

fig10.savefig(outfilepath+'source_2_some_swapped.eps')
fig9.savefig(outfilepath+'source_1_some_swapped.eps')
fig8.savefig(outfilepath+'source_2_all_swapped.eps')
fig7.savefig(outfilepath+'source_1_all_swapped.eps')
fig1.savefig(outfilepath+'orig_data.eps')
fig2.savefig(outfilepath+'orig_spur_data.eps')
fig3.savefig(outfilepath+'class_data.eps')
fig4.savefig(outfilepath+'class_pcs.eps')
fig5.savefig(outfilepath+'source_1.eps')
fig6.savefig(outfilepath+'source_2.eps')

if (show_plots):
    plt.show()
