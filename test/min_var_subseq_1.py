# Find the subsequence with minimum variance and see if it corresponds to
# something iteresting
import numpy as np

f1=100.
f2=f1*2.**(1./12.)
N=10
f=np.r_[np.arange(N)*f1+np.random.standard_normal(N),
        np.arange(N)*f2+np.random.standard_normal(N),
        np.random.uniform(1,2000.,N)]
print f
d=np.subtract.outer(f,f).flatten()
d=d[np.where(abs(d) < 300.)]
d=d[np.where(abs(d) != 0.)]
d=np.sort(d)
print d
md=max(abs(d))
mdk=0
for k in xrange(len(d)-N+1):
    md_=abs(d[k+N-1]-d[k])
    if md_ < md:
        md=md_
        mdk=k
print d[mdk:mdk+N]
