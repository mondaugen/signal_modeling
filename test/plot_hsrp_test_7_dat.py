import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

show_plots=True
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
for i in xrange(len(datin)):
    X_plt_orig=datin[i]['X_plt_orig'][0][0]
    w0_o=X_plt_orig[:,0]+(-0.5)*X_plt_orig[:,1]*H
    w1_o=X_plt_orig[:,0]+(0.5)*X_plt_orig[:,1]*H
    N_w0_o=len(w0_o)
    h_o=h+H*np.c_[-0.5*np.ones(N_w0_o),0.5*np.ones(N_w0_o)]
    ax1.plot(h_o.T[:,:K],np.c_[w0_o[:K],w1_o[:K]].T,ls_s1)
    ax1.plot(h_o.T[:,K:],np.c_[w0_o[K:],w1_o[K:]].T,ls_s2)
#    ax1.legend()
    h+=H

h=0
for i in xrange(len(datin)):
    X_plt=datin[i]['X_plt'][0][0]
    w0=X_plt[:,0]+(-0.5)*X_plt[:,1]*H
    w1=X_plt[:,0]+(0.5)*X_plt[:,1]*H
    N_w0=len(w0)
    h_=h+H*np.c_[-0.5*np.ones(N_w0),0.5*np.ones(N_w0)]
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
    for r in np.c_[h_e,w0_e,w1_e,clr_e]:
        # Source 1 'b'
        if (r[4]=='b'):
            ls_=ls_s1
        # Source 2 'c'
        elif (r[4]=='c'):
            ls_=ls_s2
        else:
            ls_=ls_bg
        ax3.plot(r[0:2],r[2:4],ls_)
    h+=H

h=0
for i in xrange(len(datin)):
    A_pca=datin[i]['A_pca'][0][0]
    clr_e=datin[i]['clr_'][0][0]
    a0=A_pca[:,0]
    N_a0=len(a0)
    h_e=h+H*np.c_[np.ones(N_a0)]
    for r in np.c_[h_e,a0,clr_e]:
        # Source 1 'b'
        if (r[2]=='b'):
            ls_s=ps_s1_s
            ls_c=ps_s1_c
            lw_=0
        # Source 2 'c'
        elif (r[2]=='c'):
            ls_s=ps_s2_s
            ls_c=ps_s2_c
            lw_=0
        else:
            ls_s=ps_bg_s
            ls_c=ps_bg_c
            lw_=0
        ax4.scatter(r[0],r[1],marker=ls_s,c=ls_c,lw=lw_)
    h+=H

if (show_plots):
    plt.show()
