% Matrix of all possible partitionings
P=16;
K=2^P;
C=zeros(K,1);
f1=rand(1)*(1:P/2);
f2=rand(1)*(1:P/2);
f=[f1 f2];
f=sort(f);
for k=0:(K-1)
    A=dec2bin(k,P)-48;
    idx1=find(A);
    idx2=find(~A);
    d1=diff(f(idx1));
    if (length(d1)==0)
        d1=0;
    end
    d2=diff(f(idx2));
    if (length(d2)==0)
        d2=0;
    end
    C(k+1)=var(d1)+var(d2);
end
display('Estimated');
[minc,minci]=min(C)
A=dec2bin(minci-1,P)-48
f(find(A))
f(find(~A))
display('True');
f1
f2
