# Plot the partial tracks computed using test/partial_estimation_plot_3_3b.py
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sigmod as sm

with open('tmp/ac_gtr_a3_op_fs3_pizz_sr16k.pardat','r') as f:
    p_info=pickle.load(f)

mpl.rcParams['legend.fontsize'] = 10

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
fig2 = plt.figure(2)
ax2= fig2.gca()
X=[]
for p in p_info:
    if len(p) > 0:
        tpts=[]
        fpts=[]
        apts=[]
        for i in xrange(len(p)):
            tpts.append(p[i][0])
            fpts.append(float(np.imag(p[i][1][1])))
            apts.append(float(np.real(p[i][1][0])))
        # fit line to data
        th_a=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],apts)[0]
        th_f=np.linalg.lstsq(np.c_[np.ones(len(tpts)),tpts],fpts)[0]
        ax.plot(tpts,fpts,apts,c='g')
        tpts=np.array(tpts)
        ax.plot(tpts,th_f[0]+th_f[1]*tpts,th_a[0]+th_a[1]*tpts,c='b')
#        ax2.scatter(th_f[0],th_a[1])
#        X.append([len(fpts),1./th_f[0]])
        X.append([len(fpts),1./th_f[0]])


A=sm.pca_ne(np.array(X).T,'cov')
#ax2.scatter(A[0,:],A[1,:])
X=np.array(X).T
ax2.scatter(X[0,:],X[1,:])

plt.show()
