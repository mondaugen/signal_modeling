clear;
N=10;
M=2*(N-1);
P=5;
noisev=1e-5;
f=rand(P,1);
fs_orig=f*(1:N)+randn(P,N)*sqrt(noisev);
fs=sort(fs_orig(:));
K=length(fs);
df=fs-fs';
df=abs(df(:));
[sdf,sdfi]=sort(df);
L=length(df);
C=zeros(L-M+1,1);
for(l=1:(L-M+1))
    if(sdf(l) > 0)
        C(l)=var(sdf(l:(l+M-1)));
    else
        C(l)=Inf;
    end
end
[cmins,cminis]=sort(C);
R=3*P;
results=zeros(M,R);
for p=1:R
    g=sdfi(cminis(p):(cminis(p)+M-1))(:);
    G=[mod(g-1,K)+1 floor((g-1)/K)+1];
    Gu=unique(G);
%    length([fs(Gu);zeros(M-length(Gu),1)])
    results(:,p)=[fs(Gu);zeros(M-length(Gu),1)];
end
[results;sum(results>0)]
fs_orig
