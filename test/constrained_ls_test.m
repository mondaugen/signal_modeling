% Solve constrained linear least squares problem.

% Closed form:

% Number of rows
m=4;
% Number of columns
n=3;
% Number of rows in contstraint matrix A
p=1;
% observation matrix
H=randn(m,n);
% observations
x=randn(m,1);
% constraint matrix
A=randn(p,n);
% constraint vector
b=randn(p,1);
% Unconstrained solution
th=linsolve(H,x);
% Constrained solution
th_c=linsolve(H'*H,H'*x-A'*(A*(H'*H)^(-1)*A')^(-1)*(A*th-b));
th_c

% Solution via QR factorization
% We use the approach described in Golub and Van Loan, pp. 266-7
B=blkdiag(eye(m),zeros(p));
C=[H;A];
c=[x;b];
[m_C,n_C]=size(C);
[Q,R]=qr(C);
Q1=Q(:,1:n_C);
Q2=Q(:,n_C+1:end);
R1=R(1:n_C,:);
% In text, they say form 
%
% Q2'*B*Z=[0,S]
%
% where S is upper-diagonal
% To use the qr factorization here we do
%
% Q2'*B=[0,S]*Z' (Z is orthogonal matrix)
% B'*Q2=Z*[0;S']
%
% You can sort of see that flipping [0;S'] vertically, then horizontally will
% give a matrix [R;0], where R is upper triangular so we compute the QR
% factorization on Q2'*B flipped vertically then horizontally.
% The relationship is then
% Z=fliplr(flipud(Q_))
% [0,S]=flipud(fliplr(R_))'

[Q_,R_]=qr(fliplr(flipud(B'*Q2)));
Z=fliplr(flipud(Q_));
OS=flipud(fliplr(R_))';
S=OS(:,n_C+1:end);

Z1=Z(:,1:n_C);
Z2=Z(:,n_C+1:end);
u=linsolve(S,Q2'*c);
v=Z2*u;
y=linsolve(R1,Q1'*c-Q1'*B*v);
y
norm(th_c-y)
