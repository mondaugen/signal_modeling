import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sigmod as sm
import sklearn.mixture
import os

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='partial_classification_acgtr_xylo_sep_'

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# Load partials from file 1
with open('tmp/ac_gtr_a3_op_sr16k.pardat','r') as f:
    p_info1=pickle.load(f)
# Load partials from file 2
with open('tmp/xylo_fs4_sr16k.pardat','r') as f:
    p_info2=pickle.load(f)

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

# Space to store combined partial paramters
X=[]
# Combine parameters
p_info=p_info1+p_info2
# Space to store the starting amplitudes of partials
a_starts=[]
# Space to store the average frequency of partials
f_avg_calc=[]
for p in p_info:
    if len(p) > 0:
        ax4.scatter(np.imag(p[0][1][1]),np.real(p[0][1][0]))
        a_starts.append(np.real(p[0][1][0]))
        f_avg_calc.append(np.imag(p[0][1][1]))
ax4.set_title('Starting amplitude vs frequency')
ax4.set_xlabel('Frequency (radians/s)')
ax4.set_ylabel('Amplitude (log(A))')

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
            as_v_f_th[0]+as_v_f_th[1]*max(f_avg_calc)*asvf_th])

# plot the partials
len_pi1=0
for p in p_info1:
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
        ax.plot(tpts,th_f[0]+th_f[1]*tpts,th_a[0]+th_a[1]*tpts,c='b')
        # Plot linearly interpolated data (2D)
        ax3.plot(tpts/float(H),th_f[0]+th_f[1]*tpts,'b')
        X.append([float(len(fpts)),th_f[0]])
        len_pi1+=1

len_pi2=0
for p in p_info2:
    if (len(p) <= 0):
        continue
    _f=np.imag(p[0][1][1])
    _a=np.real(p[0][1][0])
    # (See above for explanation)
    if ((as_v_f_th[0]+as_v_f_th[1]*_f)*asvf_th<_a):
        tpts=[]
        fpts=[]
        apts=[]
        for i in xrange(len(p)):
            tpts.append(p[i][0])
            fpts.append(float(np.imag(p[i][1][1]))/(2.*np.pi)*Fs)
            apts.append(float(np.real(p[i][1][0])))
        th_a=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],apts)[0]
        th_f=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],fpts)[0]
#        ax.plot(tpts,fpts,apts,c='g')
        tpts=np.array(tpts)
        ax.plot(tpts,th_f[0]+th_f[1]*tpts,th_a[0]+th_a[1]*tpts,c='r')
        ax3.plot(tpts/512,th_f[0]+th_f[1]*tpts,'r')
        X.append([float(len(fpts)),th_f[0]])
        len_pi2+=1

ax.set_xlabel('Time (samples)')
ax.set_ylabel('Frequency (Hz)')
ax.set_zlabel('Amplitude (log(A))')
ax.set_title('Partial trajectories')

ax3.set_xlabel('Time (samples)')
ax3.set_ylabel('Frequency (Hz)')
ax3.set_title('Partial trajectories')


X=np.array(X).T
A=sm.pca_ne(np.c_[np.log(X[0,:]),X[1,:]].T,'cov')
ax2.scatter(A[0,:len_pi1],A[1,:len_pi1],c='k',lw=0,label='Guitar')
ax2.scatter(A[0,len_pi1:],A[1,len_pi1:],c='grey',lw=0,label='Xylophone')
ax2_title='Log-partial-length vs. frequency: principal components'
ax2.set_xlabel('$a_{0}$')
ax2.set_ylabel('$a_{1}$')
ax2.legend(loc='upper left')
fig2.savefig(plotoutpath+'true_memberships.eps')
with open(plotoutpath+'true_memberships.txt','w') as f:
    f.write(ax2_title+'%')

# Convolve data points with kernels
A_x=np.linspace(A[0,:].min(),A[0,:].max(),num=100)
A_y=np.linspace(A[1,:].min(),A[1,:].max(),num=100)
A_diff_x=(A[0,:].max()-A[0,:].min())
A_diff_y=(A[1,:].max()-A[1,:].min())
A_var_scl=1.2
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
print a_lma_ma,a_lma_mi

#np.exp(-((A_y[j]+4.)**2.+A_x[i]+3000.)**2.)#_krnl(np.array([A_y[j],A_x[i]]),A.T)
ax2.contour(A_X,A_Y,A_Z)
ax2.plot(A_x[a_lma_arg_c[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],
        A_y[a_lma_arg_r[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],'kx')


# Classify partials using amplitude modulation
gmm_means=np.c_[A_x[a_lma_arg_c[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],
        A_y[a_lma_arg_r[np.r_[a_lma_ma_arg,a_lma_mi_arg]]]]
gmm=sklearn.mixture.GMM(n_components=2,init_params='c')
gmm.weights_=np.array([a_lma_ma/(a_lma_ma+a_lma_mi),a_lma_mi/(a_lma_ma+a_lma_mi)])
gmm.means_=gmm_means#np.array([[-3000,-5.5],[-1000,-3]])
gmm_grps=gmm.fit_predict(A.T)
m1_idx=np.where(gmm_grps==0)[0]
m2_idx=np.where(gmm_grps==1)[0]
ax5.scatter(A[0,m1_idx],A[1,m1_idx],c='g')
ax5.scatter(A[0,m2_idx],A[1,m2_idx],c='b')
ax5.set_title('Estimated memberships')
ax5.set_xlabel('1st PC')
ax5.set_ylabel('2nd PC')

plt.show()
