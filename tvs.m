function [z_,obj,info]=tvs(x,delta)
% [x_]=tvs(x,delta)
%
% Total variation smoothing.
% find x_, that minimizes
%
% 0.5*norm(x_-x,2)+delta*norm(D*X_,1)
% 
% Where delta is a weighting factor.
N=length(x);
z=[zeros(N,1);
   ones(N,1)];
H=[eye(N) zeros(N);
   zeros(N) zeros(N)];
H=sparse(H);
q=[-x;delta*ones(N,1)];
D=diag(ones(N,1))-diag(ones(N-1,1),1);
D=sparse(D);
%D=diag(2*ones(N,1))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);
D(end,end-1)=-1;
%A_in=[ D -eye(N);
%      -D -eye(N);
%       zeros(N) -eye(N)];
A_in=[ D -eye(N);
      -D -eye(N);];
A_in=sparse(A_in);
lb=[-Inf*ones(N,1);zeros(N,1);];
ub=[Inf*ones(2*N,1)];
%A_lb=-Inf*ones(3*N,1);
%A_ub=zeros(3*N,1);
A_lb=-Inf*ones(2*N,1);
A_ub=zeros(2*N,1);
[z_,obj,info]=qp(z,H,q,[],[],lb,ub,A_lb,A_in,A_ub,struct('MaxIter',2000));
