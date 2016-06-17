import numpy as np
import matplotlib.pyplot as plt
import ptpath
import ptpath_test
from cvxopt import solvers
import os

outpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
fig1_outpath=outpath+'small_graph_ex.eps'
fig2_outpath=outpath+'small_graph_ex_greedy_paths.eps'
fig3_outpath=outpath+'small_graph_ex_lp_paths.eps'

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

nodes=[
        ((0,0),'0'),
        ((0,4),'1'),
        ((0,5.5),'2'),
        ((1,2),'3'),
        ((1,4),'4'),
        ((1,5),'5'),
        ((1,7),'6'),
        ((2,3),'7'),
        ((2,4),'8')
        ]

plt.figure(1)
plt.title('Possible graph connections')
plt.xlabel("$\\theta_{0}$")
plt.ylabel("$\\theta_{1}$")

ax1=plt.gca()

for node in nodes:
    plt.scatter(node[0][0],node[0][1],c='w',linewidths=0)
    ax1.annotate(
            node[1],
            xy=node[0],
            xycoords='data',
            textcoords='offset points',
            xytext=(0,0),
            ha='center',
            va='center')


F=[nodes[:3],nodes[3:7],nodes[7:]]    

# Space around the point
rds=0.15

for k in xrange(len(F)-1):
    for i in xrange(len(F[k])):
        for j in xrange(len(F[k+1])):
            x0=float(F[k][i][0][0])
            x1=float(F[k+1][j][0][0])
            y0=float(F[k][i][0][1])
            y1=float(F[k+1][j][0][1])
            ax1.annotate('',
                    xy=(x0,y0),
                    xytext=(x1,y1),
                    arrowprops=dict(
                        arrowstyle='-',
                        shrinkA=10,
                        shrinkB=10)) 

ax1.set_xticks(np.array(list(xrange(F[-1][0][0][0]+1)),dtype=np.int))
plt.savefig(fig1_outpath)

# Find distances using greedy method
D_grd=np.zeros((len(F[0]),len(F[1]),len(F[2])))
for i_ in xrange(len(F[0])):
    for j_ in xrange(len(F[1])):
        for k_ in xrange(len(F[2])):
            x0,y0=F[0][i_][0]
            x1,y1=F[1][j_][0]
            x2,y2=F[2][k_][0]
            d1=((x1-x0)**2.+(y1-y0)**2.)**0.5
            d2=((x2-x1)**2.+(y2-y1)**2.)**0.5
            D_grd[i_,j_,k_]=d1+d2

I_=list(xrange(len(F[0])))
J_=list(xrange(len(F[1])))
K_=list(xrange(len(F[2])))

# Find J shortest paths using greedy method
J=2
min_ijks=[]
cs_grd=[]
for j in xrange(J):
    min_cost=float('inf')
    min_ijk=(-1,-1,-1)
    for i_ in I_:
        for j_ in J_:
            for k_ in K_:
                if D_grd[i_,j_,k_] < min_cost:
                    min_cost=D_grd[i_,j_,k_]
                    min_ijk=(i_,j_,k_)
    min_ijks.append(min_ijk)
    cs_grd.append(min_cost)
    i_,j_,k_=min_ijk
    I_.remove(i_)
    J_.remove(j_)
    K_.remove(k_)

plt.figure(2)
plt.title('Two shortest paths using the greedy method')
plt.xlabel('$\\theta_{0}$')
plt.ylabel('$\\theta_{1}$')
ax2=plt.gca()

for node in nodes:
    plt.scatter(node[0][0],node[0][1],c='w',linewidths=0)
    ax2.annotate(
            node[1],
            xy=node[0],
            xycoords='data',
            textcoords='offset points',
            xytext=(0,0),
            ha='center',
            va='center')


F=[nodes[:3],nodes[3:7],nodes[7:]]    

for j in xrange(J):
    i_,j_,k_=min_ijks[j]
    x0,y0=F[0][i_][0]
    x1,y1=F[1][j_][0]
    x2,y2=F[2][k_][0]
    ax2.annotate('',
            xy=(x0,y0),
            xytext=(x1,y1),
            arrowprops=dict(
                arrowstyle='-',
                shrinkA=10,
                shrinkB=10,
                color='k')) 
    ax2.annotate('',
            xy=(x1,y1),
            xytext=(x2,y2),
            arrowprops=dict(
                arrowstyle='-',
                shrinkA=10,
                shrinkB=10,
                color='k')) 
print 'Total cost (greedy): %f' % (sum(cs_grd),)
ax2.set_xticks(np.array(list(xrange(F[-1][0][0][0]+1)),dtype=np.int))
plt.savefig(fig2_outpath)

# Find J shortest paths using LP
def _lp_cost_fun(a,b):
    x0,y0=a.value[0]
    x1,y1=b.value[0]
    return ((x1-x0)**2.+(y1-y0)**2.)**0.5

n_node=0
n_frame=0
# Build frames and graph
F_lp=[]
S_lp=dict()
for f_ in F:
    F_lp.append([])
    for k_f_ in xrange(len(f_)):
        f__=f_[k_f_]
        F_lp[-1].append(n_node)
        S_lp[n_node]=ptpath.LPNode(value=f__,frame_num=n_frame)
        n_node+=1
    n_frame+=1
for l_ in xrange(len(F)-1):
    for f in F_lp[l_]:
        S_lp[f].out_nodes=F_lp[l_+1]
    for f in F_lp[l_+1]:
        S_lp[f].in_nodes=F_lp[l_]
# Build linear program
d=ptpath.g_f_2lp(S_lp,F_lp,J,_lp_cost_fun,{'calc_mean':0,'min_mean_dev':0})
# Solve LP
sol=solvers.lp(d['c'],d['G'],d['h'],d['A'],d['b'])
# Extract paths
paths=ptpath_test.lp_sol_extract_paths(sol['x'],S_lp,F_lp)

plt.figure(3)
plt.title('Two shortest paths using the LP method')
plt.xlabel('$\\theta_{0}$')
plt.ylabel('$\\theta_{1}$')
ax3=plt.gca()

for node in nodes:
    plt.scatter(node[0][0],node[0][1],c='w',linewidths=0)
    ax3.annotate(
            node[1],
            xy=node[0],
            xycoords='data',
            textcoords='offset points',
            xytext=(0,0),
            ha='center',
            va='center')

for path in paths:
    for k_ in xrange(len(path)-1):
        x0,y0=S_lp[path[k_]].value[0]
        x1,y1=S_lp[path[k_+1]].value[0]
        ax3.annotate('',
                xy=(x0,y0),
                xytext=(x1,y1),
                arrowprops=dict(
                    arrowstyle='-',
                    shrinkA=10,
                    shrinkB=10,
                    color='k')) 
ax3.set_xticks(np.array(list(xrange(F[-1][0][0][0]+1)),dtype=np.int))
print 'Total cost (LP): %f' % (sol['primal objective'])
plt.savefig(fig3_outpath)

plt.show()
