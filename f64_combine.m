% Combine f64_files by adding them together
% output written to stdout.
args=argv();
x=[];
for n=1:nargin
    fi=fopen(args{n},'r');
    x_=fread(fi,Inf,'float64');
    fclose(fi);
    if (length(x_)>length(x))
        x=[x;zeros(length(x_)-length(x),1)];
    elseif (length(x)>length(x_))
        x_=[x_;zeros(length(x)-length(x_),1)];
    end
    x+=x_;
end
fwrite(stdout,x,'float64');
