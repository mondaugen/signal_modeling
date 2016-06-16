import numpy as np
import matplotlib.pyplot as plt
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
            #o_=y1-y0
            #a_=x1-x0
            #alph=np.arctan2(o_,a_)
            ##y0+=rds*np.sin(alph)
            ##x0+=rds*np.cos(alph)
            ##y1-=rds*np.sin(alph)
            ##x1-=rds*np.cos(alph)
            #l=plt.plot([x0,x1],[y0,y1],c='k')
            #for l_ in l:
            #    print l_.get_xdata(False)
            #    print l_.get_xdata(True)
            ax1.annotate('',
                    xy=(x0,y0),
                    xytext=(x1,y1),
                    arrowprops=dict(
                        arrowstyle='-',
                        shrinkA=10,
                        shrinkB=10)) 

ax1.set_xticks(np.array(list(xrange(F[-1][0][0][0]+1)),dtype=np.int))
#print ax1.get_xticks()
#print type(ax1.get_xticks())
plt.show()
