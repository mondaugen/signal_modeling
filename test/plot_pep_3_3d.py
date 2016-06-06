# Plot the partial tracks computed using test/partial_estimation_plot_3_3b.py
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sigmod as sm
import sklearn.mixture

with open('tmp/ac_gtr_a3_op_sr16k.pardat','r') as f:
    p_info1=pickle.load(f)
with open('tmp/xylo_fs4_sr16k.pardat','r') as f:
    p_info2=pickle.load(f)

Fs=16000

mpl.rcParams['legend.fontsize'] = 10

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
X=[]

p_info=p_info1+p_info2
a_starts=[]
f_avg_calc=[]
for p in p_info:
    if len(p) > 0:
        ax4.scatter(np.imag(p[0][1][1]),np.real(p[0][1][0]))
        a_starts.append(np.real(p[0][1][0]))
        f_avg_calc.append(np.imag(p[0][1][1]))

## Calculate line function for thresholding initial amplitude values of partials
## Adjust level of final line function
asvf_th=1.2
A_as_v_f=np.c_[np.ones(len(f_avg_calc)),f_avg_calc]
b_as_v_f=np.array(a_starts)
as_v_f_th=np.linalg.lstsq(A_as_v_f,b_as_v_f)[0]
ax4.plot([0,max(f_avg_calc)],
        [as_v_f_th[0]*asvf_th,
            as_v_f_th[0]+as_v_f_th[1]*max(f_avg_calc)*asvf_th])
#A_as_v_f=np.c_[np.ones(len(a_starts)),a_starts]
#b_as_v_f=np.array(f_avg_calc)
#as_v_f_th=np.linalg.lstsq(A_as_v_f,np.log(b_as_v_f))[0]
#ax4.plot(np.exp(as_v_f_th[0]+as_v_f_th[1]*(0-np.arange(-min(a_starts)))*asvf_th),
#    (0-np.arange(-min(a_starts))))

len_pi1=0
for p in p_info1:
    if (len(p) <= 0):
        continue
    _f=np.imag(p[0][1][1])
    _a=np.real(p[0][1][0])
#    if (np.exp(as_v_f_th[0]+as_v_f_th[1]*_a*asvf_th)<_f):
    if ((as_v_f_th[0]+as_v_f_th[1]*_f)*asvf_th<_a):
        tpts=[]
        fpts=[]
        apts=[]
        for i in xrange(len(p)):
            tpts.append(p[i][0])
            fpts.append(float(np.imag(p[i][1][1]))/(2.*np.pi)*Fs)
            apts.append(float(np.real(p[i][1][0])))
        # fit line to data
        th_a=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],apts)[0]
        th_f=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],fpts)[0]
#        ax.plot(tpts,fpts,apts,c='g')
        tpts=np.array(tpts)
        ax.plot(tpts,th_f[0]+th_f[1]*tpts,th_a[0]+th_a[1]*tpts,c='b')
        ax3.plot(tpts/512,th_f[0]+th_f[1]*tpts,'b')
        X.append([float(len(fpts)),th_f[0]])
        len_pi1+=1

len_pi2=0
for p in p_info2:
    if (len(p) <= 0):
        continue
    _f=np.imag(p[0][1][1])
    _a=np.real(p[0][1][0])
#    if (np.exp(as_v_f_th[0]+as_v_f_th[1]*_a*asvf_th)<_f):
    if ((as_v_f_th[0]+as_v_f_th[1]*_f)*asvf_th<_a):
        tpts=[]
        fpts=[]
        apts=[]
        for i in xrange(len(p)):
            tpts.append(p[i][0])
            fpts.append(float(np.imag(p[i][1][1]))/(2.*np.pi)*Fs)
            apts.append(float(np.real(p[i][1][0])))
        # fit line to data
        th_a=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],apts)[0]
        th_f=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],fpts)[0]
#        ax.plot(tpts,fpts,apts,c='g')
        tpts=np.array(tpts)
        ax.plot(tpts,th_f[0]+th_f[1]*tpts,th_a[0]+th_a[1]*tpts,c='r')
        ax3.plot(tpts/512,th_f[0]+th_f[1]*tpts,'r')
        X.append([float(len(fpts)),th_f[0]])
        len_pi2+=1


#ax2.scatter(A[0,:],A[1,:])
X=np.array(X).T
A=sm.pca_ne(np.c_[np.log(X[0,:]),X[1,:]].T,'cov')
#ax2.scatter(np.log(X[0,:]),X[1,:])
ax2.scatter(A[0,:len_pi1],A[1,:len_pi1],c='g')
ax2.scatter(A[0,len_pi1:],A[1,len_pi1:],c='b')
ax2.set_title('True memberships')
#ax2.plot([A[0,:].min(),A[0,:].max()],0.5*(A[1,:].max()+A[1,:].min())*np.ones(2))
#ax5.hist(A[1,:],2,weights=np.r_[np.ones(len_pi1)/float(len_pi1),np.ones(len_pi2)/float(len_pi2)])
#ax2.scatter(np.log(X[0,:len_pi1]),X[1,:len_pi1],c='g')
#ax2.scatter(np.log(X[0,len_pi1:]),X[1,len_pi1:],c='b')

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
print A_Z[a_lma_arg_r,a_lma_arg_c]

#np.exp(-((A_y[j]+4.)**2.+A_x[i]+3000.)**2.)#_krnl(np.array([A_y[j],A_x[i]]),A.T)
ax2.contour(A_X,A_Y,A_Z)
ax2.plot(A_x[a_lma_arg_c[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],
        A_y[a_lma_arg_r[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],'kx')


# Classify partials using amplitude modulation
gmm_means=np.c_[A_x[a_lma_arg_c[np.r_[a_lma_ma_arg,a_lma_mi_arg]]],
        A_y[a_lma_arg_r[np.r_[a_lma_ma_arg,a_lma_mi_arg]]]]
gmm=sklearn.mixture.GMM(n_components=2,init_params='c')
gmm.weights_=np.array([0.65,0.35])
gmm.means_=gmm_means#np.array([[-3000,-5.5],[-1000,-3]])
gmm_grps=gmm.fit_predict(A.T)
m1_idx=np.where(gmm_grps==0)[0]
m2_idx=np.where(gmm_grps==1)[0]
ax5.scatter(A[0,m1_idx],A[1,m1_idx],c='g')
ax5.scatter(A[0,m2_idx],A[1,m2_idx],c='b')
ax5.set_title('Estimated memberships')

plt.show()
