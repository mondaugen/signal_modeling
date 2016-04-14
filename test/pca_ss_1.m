% Attempt source separation by clustering principle components
fname='/tmp/test4.f64';
[Pxn,PX,ENBW_PX,opt]=compute_all_metrics(struct('fname',fname,'ppvp_a',0.13,'K',20));
N=length(Pxm);
As=cell(N,1);
Xs=cell(N,1);

newplot(figure(1));
newplot(figure(2));
hold on;
for n=1:N
    idx=Pxn{n}.X_r_mi(Pxn{n}.X_r_ex_i);
    if(length(idx)>0)
        X=zeros(length(idx),2);
        X(:,1)=Pxn{n}.psi_r(idx)./Pxn{n}.w_r(idx);
        X(:,2)=Pxn{n}.mu_r(idx)./abs(Pxn{n}.X_r(idx));
        A=pca_ne(X,'corr');
        As{n,1}=A;
        Xs{n,1}=X;
        figure(1);
        w0=Pxn{n}.w_r(idx)-0.5*Pxn{n}.psi_r(idx)*opt.H;
        w1=Pxn{n}.w_r(idx)+0.5*Pxn{n}.psi_r(idx)*opt.H;
        scatter(n*ones(length(idx),1),X(:,1));
        figure(2);
        scatter(n*ones(length(idx),1),Pxn{n}.w_r(idx));
        plot((n*ones(length(idx),2)+[-0.5,0.5])',[w0,w1]','k');
    end
end
hold off;

