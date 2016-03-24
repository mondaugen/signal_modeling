N=10;
M=N*2;
P=3;
f=rand(3,1);
fs_orig=f*(1:N);
fs=sort(fs_orig(:));
K=length(fs);
df=fs-fs';
df=abs(df(:));
[sdf,sdfi]=sort(df);
L=length(df);
C=zeros(L-M+1,1);
subseq=zeros(L,1);
ssi=1;
l=1;
lastval=sdf(l);
subseq(ssi)=sdfi(l);
ssi+=1;
l+=1;
r=1;
subseqs=cell();
subseq_l=zeros(L,1);
ep=1e-6;
while (l<L)
    if (abs(sdf(l)-lastval)<ep)
        subseq(ssi)=sdfi(l);
%        lastval=sdf(l);
        ssi+=1;
    else
        %  don't store subsequence if it just points to 0s
        subseq=subseq(find(subseq>0));
%        df(subseq)
        if (!all(df(subseq)==0))
            subseqs{r}=subseq;
            subseq_l(r)=length(subseq);
            r+=1;
        end
        subseq=zeros(L,1);
        ssi=1;
        lastval=sdf(l);
        subseq(ssi)=sdfi(l);
        ssi+=1;
    end
    l+=1;
end

%G=[mod(g-1,K)+1 floor((g-1)/K)+1]
%Gu=unique(G)

