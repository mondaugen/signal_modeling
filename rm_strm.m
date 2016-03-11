N_ARGS = 3;
args=argv();
if (nargin != N_ARGS)
    fprintf(stderr,...
    ["usage: " program_name() " H W L\n" ...
     "H : hop size\n" ...
     "W : window type\n" ...
     "    can be 'hanning'\n" ...
     "L : window length\n" ...
     "Reads stream of values of type \"float64\" from stdin."]);
    exit(-1);
end
H=str2num(args{1});
W=args{2};
L=str2num(args{3});
x=zeros(L,1);
% instead of saving complex values, we save amplitude and phase
a=zeros(L,1);
phi=zeros(L,1);
w=zeros(L,1);
psi=zeros(L,1);
mu=zeros(L,1);
t=zeros(L,1);
% Load in enough samples to fill half of x so that the first frame's centre is
% aligned with the first sample
[tmp,count]=fread(stdin,L/2,'float64');
tmp=tmp(:);
x((L/2+1):(L/2+count))=tmp;
if (count != L/2)
    x((L/2+1+count):L)=zeros(L/2-count,1);
end
do
    [X_,w,psi,mu,t]=rm(x,W);
    a=abs(X_);
    phi=arg(X_);
    fwrite(stdout,a,'float64');
    fwrite(stdout,phi,'float64');
    fwrite(stdout,w,'float64');
    fwrite(stdout,psi,'float64');
    fwrite(stdout,mu,'float64');
    fwrite(stdout,t,'float64');
    % shift over values
    x(1:(L-H))=x((H+1):L);
    [tmp,count]=fread(stdin,H,'float64');
    tmp=tmp(:);
    x((L-H+1):(L-H+count))=tmp;
until (count != H)
% pad with zeros until buffer full
x((L-H+1+count):L)=zeros(H-count,1);
[X_,w,psi,mu,t]=rm(x,W);
a=abs(X_);
phi=arg(X_);
fwrite(stdout,a,'float64');
fwrite(stdout,phi,'float64');
fwrite(stdout,w,'float64');
fwrite(stdout,psi,'float64');
fwrite(stdout,mu,'float64');
fwrite(stdout,t,'float64');
