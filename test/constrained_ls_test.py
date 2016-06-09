import numpy as np
import scipy as sp
import scipy.linalg
import sigmod as sm

# Solve constrained linear least squares problem.

# Closed form:

# Number of rows
m=4
# Number of columns
n=3
# Number of rows in contstraint matrix A
p=1
# observation matrix
H=np.asmatrix(np.random.standard_normal((m,n)))
# observations
x=np.asmatrix(np.random.standard_normal((m,1)))
# constraint matrix
A=np.asmatrix(np.random.standard_normal((p,n)))
# constraint vector
b=np.asmatrix(np.random.standard_normal((p,1)))
# Unconstrained solution
th=np.linalg.lstsq(H,x)[0]
# Constrained solution
c=np.linalg.inv(A*np.linalg.inv(H.T*H)*A.T)
c=c*(A*th-b)
c=A.T*c
c=H.T*x-c
#th_c=np.linalg.lstsq(H.T*H,H.T*x-A.T*np.linalg.inv(A*np.linalg.inv(H.T*H)*A.T)*(A*th-b))
th_c=np.linalg.lstsq(H.T*H,c)[0]
print th_c

# Solution via QR factorization
# We use the approach described in Golub and Van Loan, pp. 266-7
B=np.asmatrix(sp.linalg.block_diag(np.eye(m),np.zeros(p)));
C=np.bmat('H;A');
c=np.bmat('x;b');
m_C,n_C=C.shape;
Q,R=np.linalg.qr(C,'complete');
Q1=Q[:,:n_C];
Q2=Q[:,n_C:];
R1=R[:n_C,:];
# In text, they say form 
#
# Q2'*B*Z=[0,S]
#
# where S is upper-diagonal
# To use the qr factorization here we do
#
# Q2'*B=[0,S]*Z' (Z is orthogonal matrix)
# B'*Q2=Z*[0;S']
#
# You can sort of see that flipping [0;S'] vertically, then horizontally will
# give a matrix [R;0], where R is upper triangular so we compute the QR
# factorization on Q2'*B flipped vertically then horizontally.
# The relationship is then
# Z=fliplr(flipud(Q_))
# [0,S]=flipud(fliplr(R_))'
Q_,R_=np.linalg.qr(np.fliplr(np.flipud(B.T*Q2)),'complete')
Z=np.fliplr(np.flipud(Q_))
OS=np.flipud(np.fliplr(R_)).T
S=OS[:,n_C:]
Z1=Z[:,:n_C];
Z2=Z[:,n_C:];
u=np.asmatrix(np.linalg.solve(S,Q2.T*c));
v=Z2*u;
y=np.linalg.solve(R1,Q1.T*c-Q1.T*B*v);
print y
print np.linalg.norm(th_c-y)
print sm.lstsq_c(H,x,A,b)
