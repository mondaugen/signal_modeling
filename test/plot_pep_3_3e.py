import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pickle
import sigmod as sm
import sklearn.mixture
import os
import matplotlib.ticker as ticker

# Note: some plots produced as pdf instead of eps to preserve transparency.
# use pdftops to convert to eps. This seems to preserve transparency whereas
# direct rendering to eps from matplotlib does not.

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='partial_classification_acgtr_xylo_'

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# translucent grey
alphagrey=colors.ColorConverter().to_rgba('grey')[:3]+(0.5,)

# Load partials from file 1
fname_ptls_out='tmp/xylo_fs4_ac_gtr_a3_sr16k.ptls'
fname_pardata_in='tmp/xylo_fs4_ac_gtr_a3_sr16k.pardata'
with open(fname_pardata_in,'r') as f:
    p_info=pickle.load(f)

# Sample rate
Fs=16000
# Hop size between analysis frames
H=512

# Plotting stuff
#mpl.rcParams['legend.fontsize'] = 10

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
fig2 = plt.figure(2)
ax2= fig2.gca()
fig3 = plt.figure(3)
ax3= fig3.gca()
fig4 = plt.figure(4)
ax4= fig4.gca()
fig5 = plt.figure(5)
ax5  = fig5.gca()
fig6 = plt.figure(6)
ax6  = fig6.gca()
fig7 = plt.figure(7)
ax7  = fig7.gca()

# Space to store combined partial paramters
X=[]
# Space to store the starting amplitudes of partials
a_starts=[]
# Space to store the average frequency of partials
f_avg_calc=[]
for p in p_info:
    if len(p) > 0:
        ax4.scatter(np.imag(p[0][1][1]),np.real(p[0][1][0]),c='k')
        a_starts.append(np.real(p[0][1][0]))
        f_avg_calc.append(np.imag(p[0][1][1]))
ax4.set_title('Starting amplitude vs frequency and thresholding boundary')
ax4.set_xlabel('Frequency (radians/s)')
ax4.set_ylabel('Log-Amplitude')

# Calculate line function for thresholding initial amplitude values of partials
# If starting values of amplitude under this value, these partials are
# discarded. This hopefully filters out body resonances and sidelobes.
# Adjust level of final line function
asvf_th=1.2
A_as_v_f=np.c_[np.ones(len(f_avg_calc)),f_avg_calc]
b_as_v_f=np.array(a_starts)
as_v_f_th=np.linalg.lstsq(A_as_v_f,b_as_v_f)[0]
ax4.plot([0,max(f_avg_calc)],
        [as_v_f_th[0]*asvf_th,
            as_v_f_th[0]+as_v_f_th[1]*max(f_avg_calc)*asvf_th],c='k')
fig4.savefig(plotoutpath+'a_vs_f_thresh.eps')


# plot the partials
# storage for partial data
ptls=[]
len_pi=0
for p in p_info:
    if (len(p) <= 0):
        continue
    _f=np.imag(p[0][1][1])
    _a=np.real(p[0][1][0])
    if ((as_v_f_th[0]+as_v_f_th[1]*_f)*asvf_th<_a):
        # only plot if starting amplitude above threshold
        tpts=[]
        fpts=[]
        apts=[]
        for i in xrange(len(p)):
            tpts.append(p[i][0])
            fpts.append(float(np.imag(p[i][1][1]))/(2.*np.pi)*Fs)
            apts.append(float(np.real(p[i][1][0])))
        # fit line to partial trajectory data
        th_a=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],apts)[0]
        th_f=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],fpts)[0]
        # Plot true data (3D)
#        ax.plot(tpts,fpts,apts,c='g')
        tpts=np.array(tpts)
        # Plot linearly interpolated data (3D)
        ax.plot(tpts,th_f[0]+th_f[1]*tpts,th_a[0]+th_a[1]*tpts,c=alphagrey,lw=1.)
        # Plot linearly interpolated data (2D)
        ax3.plot(tpts/float(H),th_f[0]+th_f[1]*tpts,'b')
#        X.append([np.log(float(len(fpts))),th_f[0],th_a[0]])
        X.append([np.log(float(len(fpts)))**2.,th_f[0]])
        ptls.append([tpts,fpts,apts,th_f,th_a,p])
        len_pi+=1


ax.set_xlabel('Time ($\\times 10^{4}$ samples)',linespacing=4)
ax.set_ylabel('Frequency (KHz)',linespacing=4)
ax.set_zlabel('Log-Amplitude',linespacing=4)
ax.set_title('Partial trajectories')
# Azimuth 54, elevation 19
ax.view_init(19,54)
# Scale labels
samp_scale=1.e4
def _samp_scale_func(x,pos):
    return '{:1.1f}'.format(x/samp_scale)
freq_scale=1.e3
def _freq_scale_func(x,pos):
    return '{:1.1f}'.format(x/freq_scale)
ticks_spec_samp=ticker.FuncFormatter(_samp_scale_func)
ticks_spec_freq=ticker.FuncFormatter(_freq_scale_func)
ax.xaxis.set_major_formatter(ticks_spec_samp)
ax.yaxis.set_major_formatter(ticks_spec_freq)
fig1.savefig(plotoutpath+'partial_trajectories.pdf')

ax3.set_xlabel('Time (samples)')
ax3.set_ylabel('Frequency (Hz)')
ax3.set_title('Partial trajectories')

X=np.array(X).T
A=sm.pca_ne(X,'cov')
ax2.scatter(A[0,:],A[1,:],c='k')
ax2.set_title('Unknown memberships')
ax2.set_xlabel('1st PC')
ax2.set_ylabel('2nd PC')

# Convolve data points with kernels
A_x=np.linspace(A[0,:].min(),A[0,:].max(),num=100)
A_y=np.linspace(A[1,:].min(),A[1,:].max(),num=100)
A_diff_x=(A[0,:].max()-A[0,:].min())
A_diff_y=(A[1,:].max()-A[1,:].min())
A_var_scl=0.2
#A_var_scl=1.1
A_var_x=(A_diff_y/(A_diff_x+A_diff_y))*A_var_scl
A_var_y=(A_diff_x/(A_diff_x+A_diff_y))*A_var_scl
A_X,A_Y=np.meshgrid(A_x,A_y)
A_Z=np.zeros((A_x.size,A_y.size)).astype('d')
for i in xrange(100):
    for j in xrange(100):
        for A_ in A.T:
            A_Z[j,i]+=np.exp(-((A_var_y*(A_y[j]-A_[1]))**2.
                +(A_var_x*(A_x[i]-A_[0]))**2.))

# Find local maxima on smoothed distribution
A_lma=(np.hstack((np.full(A_Z.shape[0],False)[:,np.newaxis],
    A_Z[:,1:]>A_Z[:,:-1])).astype('bool') &
    np.hstack((A_Z[:,:-1]>A_Z[:,1:],
        np.full(A_Z.shape[0],False)[:,np.newaxis])).astype('bool'))
A_lma&=(np.vstack((np.full(A_Z.shape[1],False),
    A_Z[1:,:]>A_Z[:-1,:])).astype('bool') &
    np.vstack((A_Z[:-1,:]>A_Z[1:,:],
        np.full(A_Z.shape[1],False))).astype('bool'))
a_lma_arg_r,a_lma_arg_c=np.where(A_lma)

# Find 2 extreme maxima
a_lma_ma_arg=A_Z[a_lma_arg_r,a_lma_arg_c].argmax()
a_lma_mi_arg=A_Z[a_lma_arg_r,a_lma_arg_c].argmin()
a_lma_ma=A_Z[a_lma_arg_r,a_lma_arg_c][a_lma_ma_arg]
a_lma_mi=A_Z[a_lma_arg_r,a_lma_arg_c][a_lma_mi_arg]

print A_Z[a_lma_arg_r,a_lma_arg_c]
print a_lma_ma,a_lma_mi

#np.exp(-((A_y[j]+4.)**2.+A_x[i]+3000.)**2.)#_krnl(np.array([A_y[j],A_x[i]]),A.T)
ax2.contour(A_X,A_Y,A_Z,cmap='Greys')
ax2.plot(A_x[a_lma_arg_c[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],
        A_y[a_lma_arg_r[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],'kx')


# Classify partials using amplitude modulation
# GMM weight adjustment parameter
# Maximum weight multiplied by gmm_w_th
gmm_w_th=1.1
gmm_means=np.c_[A_x[a_lma_arg_c[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],
        A_y[a_lma_arg_r[np.r_[a_lma_ma_arg,a_lma_mi_arg]]]]
gmm=sklearn.mixture.GMM(n_components=2,init_params='c',covariance_type='full')
#gmm.covars_=np.array([[A_var_x,A_var_y],[A_var_x,A_var_y]])
gmm.weights_=np.array([a_lma_ma*gmm_w_th/(a_lma_ma*gmm_w_th+a_lma_mi),
    a_lma_mi/(a_lma_ma*gmm_w_th+a_lma_mi)])
#gmm.weights=np.array([0.5,0.5])
gmm.means_=gmm_means
gmm_grps=gmm.fit_predict(A[:2,:].T)
print gmm.covars_


m1_idx=np.where(gmm_grps==0)[0]
m2_idx=np.where(gmm_grps==1)[0]
ax5.scatter(A[0,m1_idx],A[1,m1_idx],c='k',lw=0,label='Source 1')
ax5.scatter(A[0,m2_idx],A[1,m2_idx],c='grey',lw=0,label='Source 2')
ax5.set_title('Estimated memberships')
ax5.set_xlabel('$a_{0}$')
ax5.set_ylabel('$a_{1}$')
ax5.plot(A_x[a_lma_arg_c[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],
        A_y[a_lma_arg_r[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],'kx',
        label='Local maxima')
ax5.legend(loc='upper left')
ax5.contour(A_X,A_Y,A_Z,cmap='Greys')
fig5.savefig(plotoutpath+'estimated_memberships.eps')

# Dictionary to write partial sets to
ptls_out=dict()
for grp in gmm_grps:
    ptls_out[grp]=[]


# Plot source separated partials
for ptl,grp in zip(ptls,gmm_grps):
    if (grp==0):
        ax_=ax6
    else:
        ax_=ax7
    tpts,fpts,apts,th_f,th_a,p=ptl
    ptls_out[grp].append(p)
    ax_.plot(tpts/float(H),th_f[0]+th_f[1]*tpts,'b')

with open(fname_ptls_out,'w') as f:
    pickle.dump(ptls_out,f)

plt.show()
