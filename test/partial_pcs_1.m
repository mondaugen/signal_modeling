T=1.;
h=0.1;
f1=100.;
f2=f1*2.^(1/12);
f_fm1=3.;
f_fm2=3.1;
A_fm1=1;
A_fm2=1.01;
phi_fm1=1.;
phi_fm2=1.01;
t=0:h:T;
N=length(t);
P=20;
p=1:P;
F1=p'*f1 + p'*A_fm1*sin(2.*pi*f_fm1*t+phi_fm1);
F2=p'*f2 + p'*A_fm2*sin(2.*pi*f_fm2*t+phi_fm2);
newplot(figure(1))
figure(1)
hold on
h=plot(t,F1,'b');
set(h,'DisplayName','Source 1')
h=plot(t,F2,'g');
set(h,'DisplayName','Source 2')
hold off
title('Theoretical partial paths')
xlabel('Frequency in Hz')
ylabel('Time in seconds')
print('/tmp/partial_paths.pdf')
%legend()
dF1=2.*pi*p'*f_fm1*A_fm1*cos(2.*pi*f_fm1*t+phi_fm1);
dF2=2.*pi*p'*f_fm2*A_fm2*cos(2.*pi*f_fm2*t+phi_fm2);
rm1=dF1./F1+randn(P,N)*.01;
rm2=dF2./F2+randn(P,N)*.01;
newplot(figure(2))
figure(2)
hold on
plot(t,rm1,'b');
plot(t,rm2,'g');
hold off
title('FM/F ratio')
xlabel('Ratio')
ylabel('Time in seconds')
print('/tmp/fm_f_ratio_partial_paths.pdf')
%legend('Source 1','Source 2')
X=[rm1;rm2];
A=pca_ne(X);
newplot(figure(3))
figure(3)
hold on
h=scatter(A(1:P,1),A(1:P,2),'b');
set(h,'linewidth',1)
h=scatter(A(P+1:end,1),A(P+1:end,2),'g');
set(h,'linewidth',1)
hold off
title('Principal Components')
xlabel('PC 1')
ylabel('PC 2')
print('/tmp/pc_partial_paths.pdf')
%h=legend('Source 1','Source 2')
%set(h,'linewidth',1)
