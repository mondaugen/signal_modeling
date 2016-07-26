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
import matplotlib.ticker as ticker
import matplotlib.lines as mlines


plt.rc('text',usetex=True)
plt.rc('font',family='serif')
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

#show_plots=False
show_plots=False
save_figs=True
infilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'

# specify paths through stdin

if (len(sys.argv) != 2):
    raise Exception('Specify number of columns in plot.')
N_cols=int(sys.argv[1])


#infilepath+=sys.argv[1]
infilepaths=[line.strip() for line in sys.stdin]

N_plts=len(infilepaths)
if (N_plts % N_cols) != 0:
    raise Exception('N_cols must be factor of N_plts.')

fname_idx0=0
fname_idx1=0

# Original data
fig1,ax1=plt.subplots(N_plts/N_cols,N_cols,num=1,sharex='col',sharey='row')

# Original + spurious
fig2,ax2=plt.subplots(N_plts/N_cols,N_cols,num=2,sharex='col',sharey='row')

# Original + spurious, estimated membership
fig3,ax3=plt.subplots(N_plts/N_cols,N_cols,num=3,sharex='col',sharey='row')

# 1st PC
fig4,ax4=plt.subplots(N_plts/N_cols,N_cols,num=4,sharex='col',sharey='row')

# 1st estimated source only
fig5,ax5=plt.subplots(N_plts/N_cols,N_cols,num=5,sharex='col',sharey='row')

# 2nd estimated source only
fig6,ax6=plt.subplots(N_plts/N_cols,N_cols,num=6,sharex='col',sharey='row')

# 1st estimated source only, best path connection
fig7,ax7=plt.subplots(N_plts/N_cols,N_cols,num=7,sharex='col',sharey='row')

# 2nd estimated source only, best path connection
fig8,ax8=plt.subplots(N_plts/N_cols,N_cols,num=8,sharex='col',sharey='row')

# 1st estimated source only, J best path connection
fig9,ax9=plt.subplots(N_plts/N_cols,N_cols,num=9,sharex='col',sharey='row')

# 2nd estimated source only, J best path connection
fig10,ax10=plt.subplots(N_plts/N_cols,N_cols,num=10,sharex='col',sharey='row')

# 1st estimated source only, J best path connection
fig11,ax11=plt.subplots(N_plts/N_cols,N_cols,num=11,sharex='col',sharey='row')

# 2nd estimated source only, J best path connection
fig12,ax12=plt.subplots(N_plts/N_cols,N_cols,num=12,sharex='col',sharey='row')

# 1st estimated source only, spectrogram
fig13,ax13=plt.subplots(N_plts/N_cols,N_cols,num=13,sharex='col',sharey='row')

# 2nd estimated source only, spectrogram
fig14,ax14=plt.subplots(N_plts/N_cols,N_cols,num=14,sharex='col',sharey='row')

# 1st estimated source only, real part
fig15,ax15=plt.subplots(N_plts/N_cols,N_cols,num=15,sharex='col',sharey='row')

# 2nd estimated source only, real part
fig16,ax16=plt.subplots(N_plts/N_cols,N_cols,num=16,sharex='col',sharey='row')

# 1st estimated source only, amplitude function
fig17,ax17=plt.subplots(N_plts/N_cols,N_cols,num=17,sharex='col',sharey='row')

## 2nd estimated source only, amplitude function
fig18,ax18=plt.subplots(N_plts/N_cols,N_cols,num=18,sharex='col',sharey='row')
#ax18=fig18.gca()
# The true amplitude modulation parameter
fig19,ax19=plt.subplots(N_plts/N_cols,N_cols,num=19,sharex='col',sharey='row')

# The true frequency modulation / frequency parameter
fig20,ax20=plt.subplots(N_plts/N_cols,N_cols,num=20,sharex='col',sharey='row')

# Frequency smoothed sources
# Amplitude smoothed sources
#fig21, (ax21,ax22) =plt.subplots(2,1,sharex=True,num=21)
fig21,ax21_ax22=plt.subplots(2*N_plts/N_cols,N_cols,sharex='col',sharey='row',num=21)

# Hop size
H=256
Fs=16000
K=20

# Line styles
# source 1
ls_s1='-'
# source 2
ls_s2='-'
lc_s1='k'
lc_s2='grey'
# background
ls_bg='-'
lc_bg='lightgrey'
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
lgd_ax11=[[],[]]
lgd_ax12=[[],[]]
lgd_ax13=[[],[]]
lgd_ax14=[[],[]]
ax17_lgd=[None,None]
ax19_lgd=[None,None]
ax20_lgd=[None,None]
ax21_lgd=[None,None,None,None]
ax22_lgd=[None,None]

max_fm=float('-inf')
min_fm=float('inf')
max_am=float('-inf')
min_am=float('inf')

infotext=''

for fname_idx_ in xrange(N_plts):

    infilepath=infilepaths[fname_idx_]
    if not infilepath.endswith('.dat'):
        raise Exception('File must end with .dat')
    outfilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
    outfilepath+=infilepath[:infilepath.rfind('.dat')]+'_plot_'
    infotext+='%d %d %s\n' % ((fname_idx_/N_cols),fname_idx_%N_cols,infilepath)

    datin=sio.loadmat(os.environ['HOME']+
            '/Documents/development/masters_thesis/reports/plots/'+
            infilepath)['datout'][0]
    # Theoretical signal
    x_s1_th=np.zeros((len(datin)*H,),dtype='complex_')
    # Theoretical signal
    x_s2_th=np.zeros((len(datin)*H,),dtype='complex_')
    # Theoretical amplitude functions
    a_s1_th=np.zeros((K,len(datin)*H),dtype='float64')
    # Theoretical amplitude functions
    a_s2_th=np.zeros((K,len(datin)*H),dtype='float64')

    h=0
    for i in xrange(len(datin)):
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
        lgd_ax1[0]=ax1[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_s1,c=lc_s1)[0]
        lgd_ax1[1]=ax1[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_s2,c=lc_s2)[0]
        h+=H
    
    h=0
    for i in xrange(len(datin)):
        X_plt=datin[i]['X_plt'][0][0]
        w0=X_plt[:,0]+(-0.5)*X_plt[:,1]*H
        w1=X_plt[:,0]+(0.5)*X_plt[:,1]*H
        N_w0=len(w0)
        h_=h+H*np.c_[-0.5*np.ones(N_w0),0.5*np.ones(N_w0)]
        h_/=float(Fs)
        ax2[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_.T,np.c_[w0,w1].T,'k')
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
        for r,clr__ in zip(np.c_[h_e,w0_e,w1_e],clr_e):
            # Source 1 'b'
            if (clr__=='b'):
                ls_=ls_s1
                lc_=lc_s1
                lgd_ax3_idx=0
            # Source 2 'c'
            elif (clr__=='c'):
                ls_=ls_s2
                lc_=lc_s2
                lgd_ax3_idx=1
            else:
                ls_=ls_bg
                lc_=lc_bg
                lgd_ax3_idx=2
            lgd_ax3[lgd_ax3_idx]=ax3[(fname_idx_/N_cols),fname_idx_%N_cols].plot(r[0:2],r[2:4],ls_,c=lc_)[0]
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
        lgd_ax4[3]=ax4[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h_e[0],m0[0],marker=ls_s,c=ls_c,s=80)
        ls_s=ps_m2_s
        ls_c=ps_m2_c
        lgd_ax4[4]=ax4[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h_e[0],m0[1],marker=ls_s,c=ls_c,s=80)
        # Plot 1 standard deviation away from initial mean estimate using initial
        # variance estimate
        ls_s=ps_v1_s
        ls_c=ps_v1_c
        lgd_ax4[5]=ax4[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h_e[0],m0[0]+np.sqrt(s0[0]),marker=ls_s,c=ls_c,s=80)
        lgd_ax4[5]=ax4[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h_e[0],m0[0]-np.sqrt(s0[0]),marker=ls_s,c=ls_c,s=80)
        ls_s=ps_v2_s
        ls_c=ps_v2_c
        lgd_ax4[6]=ax4[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h_e[0],m0[1]+np.sqrt(s0[1]),marker=ls_s,c=ls_c,s=80)
        lgd_ax4[6]=ax4[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h_e[0],m0[1]-np.sqrt(s0[1]),marker=ls_s,c=ls_c,s=80)
        for r,clr__ in zip(np.c_[h_e,a0],clr_e):
            # Source 1 'b'
            if (clr__=='b'):
                ls_s=ps_s1_s
                ls_c=ps_s1_c
                lw_=0
                lgd_ax4_idx=0
                n_s1+=1
            # Source 2 'c'
            elif (clr__=='c'):
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
            lgd_ax4[lgd_ax4_idx]=ax4[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(r[0],r[1],marker=ls_s,c=ls_c,lw=lw_)
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
        for r,clr__ in zip(np.c_[h_e,w0_e,w1_e],clr_e):
            # Source 1 'b'
            if (clr__=='b'):
                ls_='k-'#ls_s1
                lgd_ax5_idx=0
                lgd_ax5[lgd_ax5_idx]=ax5[(fname_idx_/N_cols),fname_idx_%N_cols].plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
#        lgd_ax5[1]=ax5[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
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
        for r,clr__ in zip(np.c_[h_e,w0_e,w1_e],clr_e):
            # Source 2 'c'
            if (clr__=='c'):
                ls_='k-'#ls_s2
                lgd_ax6_idx=0
                lgd_ax6[lgd_ax6_idx]=ax6[(fname_idx_/N_cols),fname_idx_%N_cols].plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
#        lgd_ax6[1]=ax6[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
        h+=H
    
    # Arrays of start and end frequencies for each frame, each column corresponds to
    # a priori source
    F_lp=list()
    S_lp=dict()
    n_node_lp=0
    for i in xrange(len(datin)):
        X_plt_est=np.c_[datin[i]['X_plt'][0][0],
                        datin[i]['lA_plt'][0][0],
                        datin[i]['X'][0][0][:,1]]
        clr_e=datin[i]['clr_'][0][0]
        x0_s1=[]
        x0_s2=[]
        cls_e=np.ndarray((len(clr_e),),dtype='float64')
        for j in xrange(len(clr_e)):
            if clr_e[j]=='b':
                cls_e[j]=1.
            elif clr_e[j]=='c':
                cls_e[j]=2.
            else:
                cls_e[j]=0.
        for r in np.c_[X_plt_est,cls_e]:
            if (r[4]==1.):
                x0_s1.append(r[0:4])
            elif (r[4]==2.):
                x0_s2.append(r[0:4])
        S_lp[n_node_lp]=ptpath.LPNode(value=np.array(x0_s1),frame_num=i)
        S_lp[n_node_lp+1]=ptpath.LPNode(value=np.array(x0_s2),frame_num=i)
        F_lp.append([n_node_lp,n_node_lp+1])
        n_node_lp+=2
    for i in xrange(len(F_lp)-1):
        for j in F_lp[i]:
            S_lp[j].out_nodes=F_lp[i+1]
        for j in F_lp[i+1]:
            S_lp[j].in_nodes=F_lp[i]

    def _lpnode_dist_fun_f(a,b):
        tmp1_=a.value[:,0]+a.value[:,1]*float(H)*(b.frame_num-a.frame_num)
        rslt1=np.abs(np.subtract.outer(tmp1_,b.value[:,0]))
        min_dim=np.array(rslt1.shape).argmin()
        min_cxns=rslt1.min(min_dim)
        n_min_cxns=len(min_cxns)
        rslt1=np.sum(min_cxns)/n_min_cxns
        return rslt1

    def _lpnode_dist_fun_a(a,b):
        tmp2_=a.value[:,2]+a.value[:,3]*float(H)*(b.frame_num-a.frame_num)
        rslt2=np.abs(np.subtract.outer(tmp2_,b.value[:,2]))
        min_dim=np.array(rslt2.shape).argmin()
        min_cxns=rslt2.min(min_dim)
        n_min_cxns=len(min_cxns)
        rslt2=np.sum(min_cxns)/n_min_cxns
        return rslt2

    # Smooth frequency path
    J_lp=2
    D_lp=ptpath.g_f_2lp(S_lp,F_lp,J_lp,
            cost_func=_lpnode_dist_fun_f,
            opt={'calc_mean':0,'min_mean_dev':0})

    sol_lp=solvers.lp(D_lp['c'],D_lp['G'],D_lp['h'],D_lp['A'],D_lp['b'])
    q_lp=ptpath_test.lp_sol_extract_paths(sol_lp['x'],S_lp,F_lp)

    # last phases
    phi_s1_th=np.zeros((K,))
    h=0
    for i,q_ in zip(xrange(len(datin)),q_lp[0]):
        q_=int(round(q_))
        X_plt_est=S_lp[q_].value
        w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
        w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
        N_w0_e=len(w0_e)
        h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
        h_e/=float(Fs)
        for r in np.c_[h_e,w0_e,w1_e]:
            ls_='k-'#ls_s1
            lgd_ax7_idx=0
            lgd_ax7[lgd_ax7_idx]=ax7[(fname_idx_/N_cols),fname_idx_%N_cols].plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0][:K,:]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        n_=np.arange(H)
        wn0_=X_plt_orig[:,0]
        wn1_=X_plt_orig[:,0]+X_plt_orig[:,1]*H
        dwn_=(wn1_-wn0_)/float(H)
        ph_=phi_s1_th[:,np.newaxis]+np.multiply.outer(wn0_,n_)+np.multiply.outer(dwn_,n_**2.)
        phi_s1_th+=wn0_*H+dwn_*H**2.
        lA_plt_orig=datin[i]['lA_plt_orig'][0][0][:K]
        mu_orig=datin[i]['X_orig'][0][0][:K,1]
        psi_w_orig=datin[i]['X_orig'][0][0][:K,0]
        la_=lA_plt_orig+np.multiply.outer(mu_orig,n_)
        x_s1_th[h:h+H]+=np.exp(1j*ph_+la_).sum(0)
        a_s1_th[:,h:h+H]=np.exp(la_)
        a_s1_th_plt=np.exp(lA_plt_orig+np.multiply.outer(mu_orig,np.arange(-H/2,H/2)))
        ax17_lgd[0]=ax17[(fname_idx_/N_cols),fname_idx_%N_cols].plot(np.multiply.outer(np.arange(h,h+H)/float(Fs),np.ones(K)),
                a_s1_th_plt.T,c='k')[0]
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
#        lgd_ax7[1]=ax7[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_bg)[0]
        lgd_ax9[0]=ax9[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T,np.c_[w0_o,w1_o].T,'k-')[0]
        ax19_lgd[0]=ax19[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h/float(Fs)*np.ones(K),mu_orig,c='k',lw=0)
        ax20_lgd[0]=ax20[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h/float(Fs)*np.ones(K),psi_w_orig,c='k',lw=0)
        if mu_orig.max() > max_am:
            max_am=mu_orig.max()
        if mu_orig.min() < min_am:
            min_am=mu_orig.min()
        if psi_w_orig.max() > max_fm:
            max_fm=psi_w_orig.max()
        if psi_w_orig.min() < min_fm:
            min_fm=psi_w_orig.min()
        h+=H

    # last phases
    phi_s2_th=np.zeros((K,))
    h=0
    for i,q_ in zip(xrange(len(datin)),q_lp[1]):
        q_=int(round(q_))
        X_plt_est=S_lp[q_].value
        w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
        w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
        N_w0_e=len(w0_e)
        h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
        h_e/=float(Fs)
        for r in np.c_[h_e,w0_e,w1_e]:
            ls_='k-'#ls_s2
            lgd_ax8_idx=0
            lgd_ax8[lgd_ax8_idx]=ax8[(fname_idx_/N_cols),fname_idx_%N_cols].plot(r[0:2],r[2:4],ls_)[0]
        X_plt_orig=datin[i]['X_plt_orig'][0][0][K:,:]
        w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
        w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
        n_=np.arange(H)
        wn0_=X_plt_orig[:,0]
        wn1_=X_plt_orig[:,0]+X_plt_orig[:,1]*H
        dwn_=(wn1_-wn0_)/float(H)
        ph_=phi_s2_th[:,np.newaxis]+np.multiply.outer(wn0_,n_)+np.multiply.outer(dwn_,n_**2.)
        phi_s2_th+=wn0_*H+dwn_*H**2.
        lA_plt_orig=datin[i]['lA_plt_orig'][0][0][K:]
        mu_orig=datin[i]['X_orig'][0][0][K:,1]
        psi_w_orig=datin[i]['X_orig'][0][0][K:,0]
        la_=lA_plt_orig+np.multiply.outer(mu_orig,n_)
        x_s2_th[h:h+H]+=np.exp(1j*ph_+la_).sum(0)
        a_s2_th[:,h:h+H]=np.exp(la_)
        a_s2_th_plt=np.exp(lA_plt_orig+np.multiply.outer(mu_orig,np.arange(-H/2,H/2)))
        ax17_lgd[1]=ax17[(fname_idx_/N_cols),fname_idx_%N_cols].plot(np.multiply.outer(np.arange(h,h+H)/float(Fs),np.ones(K)),
                a_s2_th_plt.T,c='grey')[0]
        N_w0_o=len(w0_o)
        h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
        h_o/=Fs
#        lgd_ax8[1]=ax8[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_bg)[0]
        lgd_ax10[0]=ax10[(fname_idx_/N_cols),fname_idx_%N_cols].plot(h_o.T,np.c_[w0_o,w1_o].T,'k-')[0]
        ax19_lgd[1]=ax19[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h/float(Fs)*np.ones(K),mu_orig,c='grey',lw=0)
        ax20_lgd[1]=ax20[(fname_idx_/N_cols),fname_idx_%N_cols].scatter(h/float(Fs)*np.ones(K),psi_w_orig,c='grey',lw=0)
        h+=H

    # Smooth frequency path LP solution
    t=0
    path_marker_list=['o','s']
    path_sf_pos=dict()
    for f_lp in F_lp:
        k=0
        for f in f_lp:
            path_sf_pos[f]=k
            ax21_ax22[2*(fname_idx_/N_cols),fname_idx_%N_cols].scatter(t,k,c='k',lw=0,marker=path_marker_list[k])
            ax21_lgd[k]=mlines.Line2D([],[],color='white',lw=0,
                    markerfacecolor='k',marker=path_marker_list[k])
            k+=1
        t+=1

    path_sf_clrs=['k','grey']
    for t in xrange(len(q_lp[0])-1):
        for i in xrange(2):
            k0=path_sf_pos[q_lp[i][t]]
            k1=path_sf_pos[q_lp[i][t+1]]
            ax21_lgd[i+2]=mlines.Line2D([],[],color=path_sf_clrs[i])
            ax21_ax22[2*(fname_idx_/N_cols),fname_idx_%N_cols].annotate('',
                xy=(t,k0),
                xytext=(t+1,k1),
                arrowprops=dict(
                    linewidth=1.,
                    arrowstyle='-',
                    shrinkA=3,
                    shrinkB=3,
                    color=path_sf_clrs[i])) 

    # Smooth amplitude path
    J_lp=2
    D_lp=ptpath.g_f_2lp(S_lp,F_lp,J_lp,
            cost_func=_lpnode_dist_fun_a,
            opt={'calc_mean':0,'min_mean_dev':0})

    sol_lp=solvers.lp(D_lp['c'],D_lp['G'],D_lp['h'],D_lp['A'],D_lp['b'])
    q_lp=ptpath_test.lp_sol_extract_paths(sol_lp['x'],S_lp,F_lp)

    h=0
    for i,q_ in zip(xrange(len(datin)),q_lp[0]):
        q_=int(round(q_))
        X_plt_est=S_lp[q_].value
        w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
        w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
        N_w0_e=len(w0_e)
        h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
        h_e/=float(Fs)
        for r in np.c_[h_e,w0_e,w1_e]:
            ls_='k-'#ls_s1
            lgd_ax11_idx=0
            lgd_ax11[lgd_ax7_idx]=ax11[(fname_idx_/N_cols),fname_idx_%N_cols].plot(r[0:2],r[2:4],ls_)[0]
        h+=H

    h=0
    for i,q_ in zip(xrange(len(datin)),q_lp[1]):
        q_=int(round(q_))
        X_plt_est=S_lp[q_].value
        w0_e=X_plt_est[:,0]+(-0.5)*X_plt_est[:,1]*H
        w1_e=X_plt_est[:,0]+(0.5)*X_plt_est[:,1]*H
        N_w0_e=len(w0_e)
        h_e=h+H*np.c_[-0.5*np.ones(N_w0_e),0.5*np.ones(N_w0_e)]
        h_e/=float(Fs)
        for r in np.c_[h_e,w0_e,w1_e]:
            ls_='k-'#ls_s2
            lgd_ax12_idx=0
            lgd_ax12[lgd_ax8_idx]=ax12[(fname_idx_/N_cols),fname_idx_%N_cols].plot(r[0:2],r[2:4],ls_)[0]
        h+=H

    # Smooth amplitude path LP solution
    t=0
    path_sa_pos=dict()
    for f_lp in F_lp:
        k=0
        for f in f_lp:
            path_sa_pos[f]=k
            ax21_ax22[2*(fname_idx_/N_cols)+1,fname_idx_%N_cols].scatter(
                    t,k,c='k',lw=0,marker=path_marker_list[k])
            k+=1
        t+=1

    path_sa_clrs=['k','grey']
    for t in xrange(len(q_lp[0])-1):
        for i in xrange(2):
            k0=path_sa_pos[q_lp[i][t]]
            k1=path_sa_pos[q_lp[i][t+1]]
            ax21_ax22[2*(fname_idx_/N_cols)+1,fname_idx_%N_cols].annotate('',
                xy=(t,k0),
                xytext=(t+1,k1),
                arrowprops=dict(
                    linewidth=1.,
                    arrowstyle='-',
                    shrinkA=3,
                    shrinkB=3,
                    color=path_sa_clrs[i])) 

    # Reformat time ticks on spectrogram
    sgram_tscale=(2.*np.pi)/float(Fs)
    def _time_scale_func(x,pos):
        return '{:1.1f}'.format(x*sgram_tscale)
    ticks_spec_time=ticker.FuncFormatter(_time_scale_func)
    
    ax13[(fname_idx_/N_cols),fname_idx_%N_cols].specgram(x_s1_th,Fs=2.*np.pi,cmap='Greys')
    ax13[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0.,(h-H)/(2.*np.pi))
    ax13[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0.,np.pi)
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax13[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax13[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
    ax13[(fname_idx_/N_cols),fname_idx_%N_cols].xaxis.set_major_formatter(ticks_spec_time)
#    ax13[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Spectrogram of source 1 (true)')
    
    ax14[(fname_idx_/N_cols),fname_idx_%N_cols].specgram(x_s2_th,Fs=2.*np.pi,cmap='Greys')
    ax14[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0.,(h-H)/(2.*np.pi))
    ax14[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0.,np.pi)
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax14[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax14[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
    ax14[(fname_idx_/N_cols),fname_idx_%N_cols].xaxis.set_major_formatter(ticks_spec_time)
#    ax14[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Spectrogram of source 2 (true)')
    
    ax15[(fname_idx_/N_cols),fname_idx_%N_cols].plot(np.arange(len(x_s1_th))/float(Fs),np.real(x_s1_th),c='k')
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax15[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax15[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Value (real part)')
#    ax15[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 1 (true)')
    ax15[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0.,(h-H)/float(Fs))
    
    ax16[(fname_idx_/N_cols),fname_idx_%N_cols].plot(np.arange(len(x_s2_th))/float(Fs),np.real(x_s2_th),c='k')
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax16[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax16[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Value (real part)')
#    ax16[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 2 (true)')
    ax16[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0.,(h-H)/float(Fs))
    
    #ax17[(fname_idx_/N_cols),fname_idx_%N_cols].plot(np.multiply.outer(np.ones((K,)),np.arange(a_s2_th.shape[1])).T,
    #        a_s1_th.T,ls='-',c='k')
#    ax17[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Amplitude function for each source (true)')
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax17[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax17[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Amplitude')
#    ax17[(fname_idx_/N_cols),fname_idx_%N_cols].legend(ax17_lgd,('Source 1','Source 2'))
    ax17[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    #ax18[(fname_idx_/N_cols),fname_idx_%N_cols].plot(np.multiply.outer(np.ones((K,)),np.arange(a_s2_th.shape[1])).T,
    #        a_s2_th.T,ls='-',c='k')
    
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax19[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax19[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('$\\mu_{k,p}$')
#    ax19[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Amplitude modulation (true)')
    mima_am_de=(max_am-min_am)*0.1
    ax19[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(min_am-mima_am_de,max_am+mima_am_de)
    ax19[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0.,(h-H)/float(Fs))
#    ax19[(fname_idx_/N_cols),fname_idx_%N_cols].legend(ax19_lgd,('Source 1','Source 2'))
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax20[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax20[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('$\\psi_{k,p}$')
#    ax20[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Frequency modulation (true)')
    mima_fm_de=(max_fm-min_fm)*0.1
    ax20[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(min_fm-mima_fm_de,max_fm+mima_fm_de)
    ax20[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0.,(h-H)/float(Fs))
#    ax20[(fname_idx_/N_cols),fname_idx_%N_cols].legend(ax20_lgd,('Source 1','Source 2'))
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax6[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax6[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax6[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 2 (estimated)')
    ax6[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax6[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax5[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax5[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax5[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 1 (estimated)')
    ax5[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax5[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax4[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax4[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('1st PC')
#    ax4[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Principal components and their classification')
    ax4[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0.1,0.5)
    ax4[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(-0.0008,0.0008)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax3[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax3[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax3[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Classified data-points')
    ax3[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax3[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    ax2[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax2[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax2[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax2[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax2[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Original and spurious data-points')
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax1[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax1[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax1[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Original data-points')
    ax1[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax7[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax7[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax7[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 1 (estimated) after smooth frequency path search')
    ax7[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax7[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax8[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax8[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax8[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 2 (estimated) after smooth frequency path search')
    ax8[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax8[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax11[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax11[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax11[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 1 (estimated) after smooth amplitude path search')
    ax11[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax11[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax12[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax12[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax12[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 2 (estimated) after smooth amplitude path search')
    ax12[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax12[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax9[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax9[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax9[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 1 partials (true)')
    ax9[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax9[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax10[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlabel('Time (seconds)')
    if (fname_idx_%N_cols)==0:
        ax10[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylabel('Frequency ($\\frac{\\text{rad}}{\\text{s}}$)')
#    ax10[(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Source 2 partials (true)')
    ax10[(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(0,(h-H)/float(Fs))
    ax10[(fname_idx_/N_cols),fname_idx_%N_cols].set_ylim(0,np.pi)
    
#    ax1[(fname_idx_/N_cols),fname_idx_%N_cols].legend(lgd_ax1,('Source 1','Source 2'))
#    ax3[(fname_idx_/N_cols),fname_idx_%N_cols].legend(lgd_ax3,('Source 1','Source 2','Spurious'))
#    ax4[(fname_idx_/N_cols),fname_idx_%N_cols].legend(lgd_ax4,('Source 1','Source 2','Spurious','$\\mu_{0}^{0}$',
#                        '$\\mu_{1}^{0}$','$\\pm \\sigma_{0}^{0}$',
#                        '$\\pm \\sigma_{1}^{0}$'))
    
    ax21_ax22[2*(fname_idx_/N_cols),fname_idx_%N_cols].get_yaxis().set_visible(False)
    ax21_ax22[2*(fname_idx_/N_cols)+1,fname_idx_%N_cols].get_yaxis().set_visible(False)
    ax21_ax22[2*(fname_idx_/N_cols),fname_idx_%N_cols].set_xlim(-0.5,len(q_lp[0])-1+0.5)
    ax21_ax22[2*(fname_idx_/N_cols)+1,fname_idx_%N_cols].set_xlim(-0.5,len(q_lp[0])-1+0.5)
#    ax21_ax22[2*(fname_idx_/N_cols),fname_idx_%N_cols].set_title('Smoothed frequency paths')
#    ax21_ax22[2*(fname_idx_/N_cols)+1,fname_idx_%N_cols].set_title('Smoothed amplitude paths')
    if (fname_idx_/N_cols) == (N_plts/N_cols-1):
        ax21_ax22[2*(fname_idx_/N_cols)+1,fname_idx_%N_cols].set_xlabel('Frame number $h$')
    #ax21[(fname_idx_/N_cols),fname_idx_%N_cols].legend(ax21_lgd,
    #    ('Original Source 1',
    #        'Original Source 2',
    #        'Source 1',
    #        'Source 2'))
    #ax22[(fname_idx_/N_cols),fname_idx_%N_cols].legend(ax22_lgd,('Source 1','Source 2'))

outfilepath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/hsrp_test_7_multi_'

if (save_figs):
    fig1.savefig(outfilepath+ 'orig_data.eps')
    fig2.savefig(outfilepath+ 'orig_spur_data.eps')
    fig3.savefig(outfilepath+ 'class_data.eps')
    fig4.savefig(outfilepath+ 'class_pcs.eps')
    fig5.savefig(outfilepath+ 'source_1_est.eps')
    fig6.savefig(outfilepath+ 'source_2_est.eps')
    fig7.savefig(outfilepath+ 'source_1_smooth_freq.eps')
    fig8.savefig(outfilepath+ 'source_2_smooth_freq.eps')
    fig9.savefig(outfilepath+ 'source_1_true.eps')
    fig10.savefig(outfilepath+'source_2_true.eps')
    fig11.savefig(outfilepath+'source_1_smooth_amp.eps')
    fig12.savefig(outfilepath+'source_2_smooth_amp.eps')
    fig13.savefig(outfilepath+'source_1_spec.eps')
    fig14.savefig(outfilepath+'source_2_spec.eps')
    fig15.savefig(outfilepath+'source_1_tdrp.eps')
    fig16.savefig(outfilepath+'source_2_tdrp.eps')
    fig17.savefig(outfilepath+'af.eps')
    fig19.savefig(outfilepath+'mu.eps')
    fig20.savefig(outfilepath+'psi.eps')
    fig21.savefig(outfilepath+'smooth_freq_amp_sol.eps')

with open(outfilepath+'_plots.INFO','w') as f:
    f.write(infotext)

if (show_plots):
    plt.show()
