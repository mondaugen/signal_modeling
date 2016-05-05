import matplotlib.pyplot as plt
import random
x=[random.random() for _ in xrange(1000)]
y=[random.random() for _ in xrange(1000)]
f1=plt.figure(1,figsize=(7,4.5))
plt.scatter(x,y)
plt.xlim((0,1))
plt.ylim((0,1))
plt.savefig('/tmp/fig1.eps',bbox_inches='tight')
f2=plt.figure(2,figsize=(10,6.5))
plt.scatter(x,y)
plt.ylim((0,1))
plt.savefig('/tmp/fig2.eps',bbox_inches='tight')
