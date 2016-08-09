import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import os

# Illustrates ad hoc cluster paritioning with a graph

show_plot=False

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='ad_hoc_cluster_part_ex.eps'

mu0=0.
mu1=5.
mu2=10.
sig0=1.
sig1=1.
sig2=1.
sig_kern=1.**0.5
N_dp=100
x=np.r_[np.random.standard_normal(N_dp)*sig0+mu0,
        np.random.standard_normal(N_dp)*sig1+mu1,
        np.random.standard_normal(N_dp)*sig2+mu2]

x0=x.min()
x1=x.max()
xs=np.linspace(x0,x1,500)
y=np.sum(np.exp(-(np.subtract.outer(x,xs)/sig_kern)**2.),0)
y/=np.sum(y)
y_ma,y_mai=sm.lextrem(y,comp='max')
y_mi,y_mii=sm.lextrem(y,comp='min')
#plt.plot(xs,y,c='k')
plt.fill_between(xs[:y_mii[0]],0.,y[:y_mii[0]],color='lightgrey',label='Partition 1')
plt.fill_between(xs[y_mii[0]:y_mii[1]],0.,y[y_mii[0]:y_mii[1]],color='grey',label='Partition 2')
plt.fill_between(xs[y_mii[1]:],0.,y[y_mii[1]:],color='dimgrey',label='Partition 3')
#plt.scatter(xs[y_mai],y_ma,c='k',lw=0.,marker='o')
l4=plt.scatter(xs[y_mii],y_mi,c='k',lw=0.,marker='o')
plt.xlim(x0,x1)
plt.ylim(0.,np.max(y_ma)*1.1)
plt.xlabel('$\\theta$')
plt.ylabel('$\\mathrm{p}(\\theta)$')
#plt.title('Partitioning example')
with open(plotoutpath[:plotoutpath.rfind('.eps')]+'_plot_title.txt','w') as f:
    f.write('Paritioning example'+'%')

# build legend
l1=plt.Rectangle((0,0),10,10,color='lightgrey',lw=0)
l2=plt.Rectangle((0,0),10,10,color='grey',lw=0)
l3=plt.Rectangle((0,0),10,10,color='dimgrey',lw=0)
plt.legend((l1,l2,l3,l4),('Partition 1','Partition 2','Partition 3',
    'Local minima'),loc='best')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off')

plt.savefig(plotoutpath)

if (show_plot):
    plt.show()
