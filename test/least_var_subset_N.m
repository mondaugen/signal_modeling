% Size of set
P=12;
% Size of subset
N=4;
K=2^P;
C=zeros(K,1);
f=[randn(P/2,1)*rand(1);randn(P/2,1)*rand(1)];
f=sort(f);
for k=0:(K-1)
    A=dec2bin(k,P)-48;
    if (sum(A)==N)
        idx1=find(A);
        d1=f(idx1);
        C(k+1)=var(d1);
    else
        C(k+1)=Inf;
    end
end
f
[cmin,cmini]=min(C);
cmin
x_e=dec2bin(cmini-1,P)-48
x0=dec2bin(2^N-1,P)-48;
x0=x0(:);
H=diag(f.^2)-f*f'/N;
[x,obj,info]=qp(x0,H,zeros(P,1),ones(1,P),N,zeros(P,1),ones(P,1));
x'
