function f=heComputeEigs(pr,beta,D,NNbar,ntot,S)

Deff=D.*repmat(S,1,ntot)./repmat(NNbar',ntot,1);
onesn=ones(ntot,1);

F=zeros(3*ntot,3*ntot);
%F=zeros(4*ntot,4*ntot);
F(1:ntot,ntot+1:end)=[pr.red*Deff,Deff];
%F(1:ntot,ntot+1:end)=[pr.red*Deff,repmat(Deff,1,2)];

vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       (pr.g2+pr.h).*onesn];
%vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       pr.g2.*onesn;       (pr.gX+pr.h).*onesn];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=    diag(-(1-pr.p3) .*(1-pr.p1) .*pr.sigma  .*onesn);
V(2*ntot+1:3*ntot,1:ntot)=  diag(-(1-pr.p4) .*pr.p1     .*pr.sigma  .*onesn);
%V(2*ntot+1:3*ntot,1:ntot)=  diag(-(1-pr.p4) .*pr.p1     .*(1-pr.p2) .*pr.sigma  .*onesn);
%V(3*ntot+1:4*ntot,1:ntot)=  diag(-(1-pr.p5) .*pr.p1     .*pr.p2     .*pr.sigma  .*onesn);

GD=F/V;

%vvec=kron([pr.sigma;pr.g1;pr.g2;pr.gX],ones(ntot,1));
%vvec(end-ntot+1:end)=vvec(end-ntot+1:end)+pr.h;

f=eigs(GD,1);%largest in magnitude (+/-)
f=max(f);
f=beta*f;%the eigenvalues of beta*F/V are beta times the eigenvalues of F/V

end