# Find the subsequence with minimum variance and see if it corresponds to
# something iteresting
import numpy as np

N=10
f1=np.arange(1,N+1)*100.+np.random.standard_normal(N)
print f1
f2=np.arange(1,N+1)*100.*2.**(1./12.)+np.random.standard_normal(N)
print f2
f=np.r_[f1,
        f2,
        np.random.uniform(1,2000.,size=N)]
L=len(f)
#index i of d corresponds to difference f[i%L] - f[floor(i/L)]
d_=np.subtract.outer(f,f).T
d=d_.flatten()
i=np.argsort(d)
s=d[i]
k=0
while s[k] <= 0.:
    k+=1
md=max(abs(s))
mdk=k
while k < (len(s)-N+1):
    if s[k+N-1] > 300.:
        break
    md_=s[k+N-1]-s[k]
    if md_ < md:
        md=md_
        mdk=k
    k+=1
for j in xrange(N):
    r=i[mdk+j]%L
    c=int(i[mdk+j])/L
    print '(%d,%d) = %f: %f, %f' % (r,c,d_[c,r],f[c],f[r])
