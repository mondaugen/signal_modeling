function f0sc_strm_f(N,H,W,L,pch_low,pch_high,c,Fs,P,format,datapath)
% function f0sc_strm_f(N,H,W,L,pch_low,pch_high,c,Fs,P,format,datapath)
x=zeros(L,1);
switch datapath
case '-'
    fin=stdin;
    fout=stdout;
otherwise
    fin=fopen(datapath,'r');
    if (fin < 0)
        error(sprintf('Error opening %s.',datapath));
    end
    fout=fopen([datapath '.f0sc'],'w');
    if (fout < 0)
        error(sprintf('Error opening %s.',[datapath '.f0sc']));
    end
end
fmt=0;
% scaling of cepstrum at extremes of pitch range
alpha=0.5;
sgm_c=sqrt(-(pch_high-pch_low)^2/(8*log(alpha)));
F=(0:(N-1))/N*Fs;
F=F(:);
F_c=Fs./(0:(N-1));
F_c=F_c(:);
w_s=exp(-(12*log(F_c./440)/log(2)+69-(pch_low+pch_high)/2).^2/(2*sgm_c^2));
% cepstrum window hardcoded as blackman, for now
L_c=round(L*c);
w_c=[blackman(L_c);zeros(N-L_c,1)];
switch format
case "raw"
    fmt=1;
case "text"
    fmt=2;
end
if (!fmt)
    error(sprintf('Bad format %s.',format));
end
% Load in enough samples to fill half of x so that the first frame's centre is
% aligned with the first sample
[tmp,count]=fread(fin,L/2,'float64');
tmp=tmp(:);
x((L/2+1):(L/2+count))=tmp;
if (count != L/2)
    x((L/2+1+count):L)=zeros(L/2-count,1);
end
do
    [Px,ENBW]=mper(x,W,N,L);
    % cepstrum window hardcoded as blackman, for now
    s=abs(ifft(log(Px).*w_c));
    inst_pow=sum(Px)/ENBW;
    s=s.*w_s;
    if (P>=0)
        [smax,sidx]=sort(s,'descend');
        s=Fs./(sidx-1);
        s=12*log(s/440.)/log(2)+69;
        s=sort(s(1:P));
    end
    switch fmt
    case 1
        % write instantaneous power
        fwrite(fout,inst_pow,'float64');
        % write f0 scores / pitches
        fwrite(fout,s,'float64');
    case 2
        % write instantaneous power, then f0 scores / pitches
        outdata=[inst_pow(:);s(:)];
        fprintf(fout,'%f ',outdata);
        fprintf(fout,'\n');
    end
    % shift over values
    x(1:(L-H))=x((H+1):L);
    [tmp,count]=fread(fin,H,'float64');
    tmp=tmp(:);
    x((L-H+1):(L-H+count))=tmp;
until (count != H)
% pad with zeros until buffer full
x((L-H+1+count):L)=zeros(H-count,1);
[Px,ENBW]=mper(x,W,N,L);
s=abs(ifft(log(Px).*w_c));
inst_pow=sum(Px)/ENBW;
s=s.*w_s;
if (P>=0)
    [smax,sidx]=sort(s,'descend');
    s=Fs./(sidx-1);
    s=12*log(s/440.)/log(2)+69;
    s=sort(s(1:P));
end
switch fmt
case 1
    % write instantaneous power
    fwrite(fout,inst_pow,'float64');
    % write f0 scores / pitches
    fwrite(fout,s,'float64');
case 2
    % write instantaneous power, then f0 scores / pitches
    outdata=[inst_pow(:);s(:)];
    fprintf(fout,'%f ',outdata);
    fprintf(fout,'\n');
end
switch datapath
case '-'
otherwise
    fclose(fin);
    fclose(fout);
end
