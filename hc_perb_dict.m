function [W]=hc_perb_dict(f,f0,N,alph,bet,w='hanning')
% HC_PERB_DICT
%
% Make harmonic dictionary whose atoms are representations of the fundamentals
% in f0. Each fundamental consists of N atoms, corresponding to N bands,
% overlapping by bet. The window shape is defined by w. The atoms are evaluated
% at frequencies in f.
if (nargin < 6)
    w='hanning';
end
f0=f0(:);
switch (w)
    case 'hanning'
        Wb=perb_hann(f,N,bet);
    otherwise
        error('Window requested not supported.');
end
W=zeros(length(f),length(f0)*N);
n=1;
k=1;
for k=1:length(f0);
    W(:,n:(n+N-1))=harm_comb(f,f0(k),alph).*Wb;
    n=n+N;
end
