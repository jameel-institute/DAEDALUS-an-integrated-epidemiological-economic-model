function [f,g,tvec]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Xit,tvec,plotTau,data)%,pmodIn)

%Inputs up to beta are outputs from hePrepCovid19
%Xit - column vector with proportion of each sector open at each
%intervention point. 
%tvec - vector of time points including tvev(1)=t0, tvec(2)=lockdown start,
%tvec(end)=end of simulation. 
%plotTau=1 to plot output. 
%Note: this is set to 0 in any optimisation protocol to avoid a crash due
%to rendering loads of images!

%%

lt=length(tvec);

% pmod=0.5408;%lockdown delta%from fitting
% pmodIn=0.6;%post-lockdown delta
% pr.betamod=[1,pmod,pmodIn*ones(1,lt-3)];

% pmod=0.65;%lockdown delta%from fitting
% pmodIn=0.7;%post-lockdown delta
% pr.betamod=[1,pmod,pmodIn*ones(1,lt-3)];

%%
isdual=1;
solvetype=2;
numseed=7;
phi1=1; phi2=0;
eps=0;
randic=0;
tauend=1;
randFact=0;
%%
%info1: (wave/info project)
%{
mcmc=1;
age2mats=0;
y0in=0; 
t0shift=0;
t0shift=tswitch-params(end-2)+120;
%}
%{
%W2 only:
y0in=ydata(xdata<18,:);
y0in=sum(y0in,1)';
%beforeShift=120;
%t0shift=tswitch-params(end-2)+120-beforeShift;%****
%}
%%
%randFact=.3;
demog=1;
%plotTau=0;
time=(1:tauend);
%ltime=length(time);
t0=0; tend=366;
%mu=1/80;%In ODE code
%phi1=1; phi2=0;
NN0=NNrep; NN0(NNrep==0)=1;%to stop blow up
Nages=NNbar./NN0;
alpha=1;%TSIR/sub-exp parameter
thetaGeom=0;
rate=eps;
%%
if tauend>1
    %Theta:
    Prow=[0,0,0,0,.5-eps,.5+eps];
    %Prow=[0,0,0,1-eps,eps];
    %
    %DO NOT TOUCH PARAMETERS:
    lp=length(Prow);
    Prow=Prow/sum(Prow);
    %
    %New for age structure:
    if na==4
        v=[5,15,45,16];
        %vover=[1/5,1/15,1/45,1/16];%[1/5,1/14,1/45,0];
        ageProps=repelem(1./v,v);
        ageProps=repmat(ageProps,n,1);
        %vover1=NNbar(1:n)/v(1); vover2=NNbar(n+1:2*n)/v(2); vover3=NNbar(2*n+1:3*n)/v(3); vover4=NNbar(3*n+1:end)/v(4);
        vover1=ones(n,1)/v(1); vover2=ones(n,1)/v(2); vover3=ones(n,1)/v(3); vover4=ones(n,1)/v(4); 
        %}
        if solvetype==3
            error('Code not written yet')
        else
            NNall=[repmat(NNbar(1:n),1,v(1)),repmat(NNbar(n+1:2*n),1,v(2)),repmat(NNbar(2*n+1:3*n),1,v(3)),repmat(NNbar(3*n+1:end),1,v(4))].*ageProps;
        end
    elseif na==5
        v=[5,13,32,15,21];
        ageProps=repelem(1./v,v);
        ageProps=repmat(ageProps,n,1);
        vover1=ones(n,1)/v(1); vover2=ones(n,1)/v(2); vover3=ones(n,1)/v(3); vover4=ones(n,1)/v(4); vover5=ones(n,1)/v(5);
        %}
        if solvetype==3
            error('Code not written yet')
        else
            NNall=[repmat(NNbar(1:n),1,v(1)),repmat(NNbar(n+1:2*n),1,v(2)),repmat(NNbar(2*n+1:3*n),1,v(3)),repmat(NNbar(3*n+1:4*n),1,v(4)),repmat(NNbar(4*n+1:end),1,v(5))].*ageProps;
        end
    end
    %No need to re-define NNbar here
end
%%
%Create time series of interventions:
%
D=Dout;
%No intervention:
%{
tvec=[0,tend];
Dvec=D;
%}

%% ORIGINAL:
 
%HE:
%Need to add age structure in here - only uses n, not na****
%tvec=[-180,92,122,153,183,tend];%,214,366];%245,275,306,336,tend];
lc=4;%*age
adInd=3;%Age group of adults *age
%lt=length(tvec);
lx=length(Xit)/(lt-3);%Number of sectors
XitMat=reshape(Xit,lx,lt-3);
NNvec=repmat(NNbar(1:lx),1,lt-3).*XitMat;%Assumes pre-lockdown=fully open
NNworkSum=sum(NNvec,1);
NNvec(lx+1:lx+lc,1:lt-3)=repmat(NNbar(lx+1:lx+lc),1,lt-3);
NNvec(lx+adInd,:)=sum(NNbar([1:lx,lx+adInd]))-NNworkSum;
%{
lt=length(tvec);
NNvec=repmat(NNbar,1,lt-1);%Full heterogeneity
Dvec=repmat(D,[1,1,lt-1]);
NNmat=reshape(NNbar,n,na);%Include at home****
NNworkSum=sum(NNmat(1:n-1,:),1);
%}
%Manually define populations during lockdown:
NNvec=[NNbar,[NNbar(1:lx).*data.xmin.^(1/pr.a)';NNbar(lx+1:lx+lc)],NNvec];
%NNvec(lx+adInd,2)=sum(NNbar([1:lx,lx+adInd]));
%NNfrac=NNvec(1:lx,2)./NNbar(1:lx,1);%Xit for lockdown
NNvec(lx+adInd,2)=NNvec(lx+adInd,2)+sum(NNvec(1:lx,1)-NNvec(1:lx,2));
%
Dvec=repmat(D,[1,1,lt-1]);
Dvec(:,:,2)=pr.betamod(2)*heMakeDs(NNvec(:,2),data.xmin.^(1/pr.a)',data,1);%,0); %NNfrac
%if lt>3

%To test attack rates etc:
%Dvec(:,:,1)=pr.betamod(1)*heMakeDs(NNvec(:,1),ones(lx,1),datax,0);

for i=3:lt-1
    %NNfrac=NNvec(1:end-1,i)./NNbar(1:end-1,1);
    Dvec(:,:,i)=pr.betamod(i)*heMakeDs(NNvec(:,i),XitMat(:,i-2),data,i-1);%NNfrac,Xit((i-3)*ntot+4));%Can reduce contact rates here too
    %{
    %Toy example:
    factor=(i-2)/(lt-2);%Includes making lockdown matrix
    NNmati=[factor*NNmat(1:n-1,:);NNmat(n,:)+(1-factor)*NNworkSum];
    NNvec(:,i)=reshape(NNmati,ntot,1);
    Dvec(:,:,i)=heMakeDs(NNvec(:,i),factor*eye(10));%Can reduce conbtact rates here too
    %}
end
%end

%% PANDEMIC PREPAREDNESS:

% lc=4;%*age
% adInd=3;%Age group of adults *age
% lx=length(data.G);%Number of sectors
% 
% XitMat=reshape(Xit,lx,lt-2);
% NNvec=repmat(NNbar(1:lx),1,lt-2).*XitMat;%Assumes pre-lockdown=fully open
% NNworkSum=sum(NNvec,1);
% NNvec(lx+1:lx+lc,1:lt-2)=repmat(NNbar(lx+1:lx+lc),1,lt-2);
% NNvec(lx+adInd,:)=sum(NNbar([1:lx,lx+adInd]))-NNworkSum;
% NNvec=[NNbar,NNvec];
% 
% Dvec=repmat(D,[1,1,lt-1]);

% for i=2:lt-1
%     Dvec(:,:,i)=pr.betamod(i)*heMakeDs(NNvec(:,i),XitMat(:,i-1),data,i-1);
% end

%%
if tauend>1
    %INITIAL CONDITION:
    B=zeros(ntot,lp);
    if randic==1
        if solvetype==3
            if ICdepAge==1
                b=zeros(ntot,lp);
                b(:,1)=NNbar;
                for nn=1:ntot
                    if NNbar(nn)>0
                        r1=round((lp-1)*rand);%How many numbers of years non-zero
                        r1=randsample(2:lp,r1);%Which years non-zero
                        toMove=round(randFact*rand*NNbar(nn));
                        if toMove>0 && isempty(r1)==0
                            toHere=datasample(r1,toMove);
                            v=ones(size(toHere));
                            toAdd=accumarray(toHere(:),v(:));
                            b(nn,2:length(toAdd))=toAdd(2:end);
                            b(nn,1)=b(nn,1)-sum(toAdd);
                        end
                    end 
                end
            else
                    error('This bit of code isnt written yet')
            end
        else
            if ICdepAge==1
                b=zeros(ntot,lp);
                b(:,1)=NNbar;
                for nn=1:ntot
                    if NNbar(nn)>0
                        r1=round((lp-1)*rand);%How many numbers of years non-zero
                        r1=randsample(2:lp,r1);%Which years non-zero
                        toMove=randFact*rand*NNbar(nn);
                        if toMove>0 && isempty(r1)==0
                            nnvec=rand(1,length(r1));
                            nnvec=nnvec/sum(nnvec)*toMove;
                            b(nn,r1)=nnvec;
                            b(nn,1)=b(nn,1)-toMove;
                        end
                    end 
                end
            else
                error('This bit of code needs checking')
                %Need to rewrite in terms of B
                maxYears=lp-1;
                Prows=zeros(n,maxYears);
                numNonZero=floor((maxYears+1)*rand(n,1));%0-6 years - pick whow many are non zero
                for i=1:n
                    rsi=randsample(maxYears,numNonZero(i));
                    Prows(i,rsi)=1;
                end
                findNZ=find(Prows);
                lNZ=length(findNZ);
                randNZ=rand(lNZ,1);
                Prows(findNZ)=randNZ;
                rand2=rand(n,1); repRand=repmat(rand2,1,maxYears);
                sumP=sum(Prows,2); repP=repmat(sumP,1,maxYears);
                B=Prows./repP.*repRand;
                B(isinf(B)==1)=0; B(isnan(B)==1)=0;
                NNratio=repmat(NNbar./NNrep,1,maxYears);
                B=[zeros(ntot,1),repmat(B,4,1).*NNratio];
            end
        end
    else
        b=zeros(ntot,lp);
        b(:,1)=NNbar;
    end
    S0=b(:,1);
    Z0=sum(b(:,2:end),2);%./NNbar;%*(1-mu);
else
    Z0=zeros(ntot,1);%Option - as input
end
%%
A1=zeros(n,lt);
A2=A1;
seed=10^-(numseed);%*NNprob;
seedvec=seed*ones(ntot,1); %zeros(ntot,1); seedvec(2*n+1:3*n)=seed*ones(n,1);
%seedvec=seed*ones(ntot,1);
thresh=0;%Remove from ODE solver
%%
%SIMULATE (UP TO ATTACK RATES):
for t=1:tauend%t=tau
[DEout,Rout,tvec]=heSimCovid19vax(pr,vx,beta,tvec,Dvec,n,ntot,NNvec,phi1,phi2,seedvec,NNvec(:,1),t,plotTau);%S0=NNbar (2nd last arg)
%[DEout,Rout]=heSimCovid19vax(pr,beta,tvec,Dvec,n,ntot,NNvec,phi1,phi2,seedvec,NNvec(:,1),t,plotTau);%S0=NNbar (2nd last arg)
%nu=Zsol-Z0;
A1=DEout;%(:,t)=sum(reshape(Zsol./NN0,n,na),2);%Prop immune for spatial cell (before antigenic drift)
A2=Rout;%(:,t)=sum(reshape(nu./NN0,n,na),2);%AR for spatial cell
%%
if tauend>1
    %ASSIGNING IMMUNITY:
    if solvetype==3
        badd=zeros(ntot,lp);
        b(:,1)=b(:,1)-nu;
        for ii=1:lp-1%Short loop - check repeated use of binomial
            baddii=binornd(nu,Prow(ii));%Were sus to both, now sus to s2/immune to just s1
            badd(:,ii)=badd(:,ii)+baddii;
            nu=nu-baddii;
        end
        badd(:,end)=nu;
        %b=b+badd;
    else
        b(:,1)=b(:,1)-nu;
        %b=b+repmat(nu,1,lp).*Prow;
    end
    %%
    %DRIFT:
    b=[b(:,1)+b(:,2),b(:,3:end),zeros(ntot,1)];
    %%
    if solvetype==3
        b=b+badd;
    else
        b=b+repmat(nu,1,lp).*repmat(Prow,ntot,1);
    end
    %%
    %AGEING:
    if demog==1
        if solvetype==3
            b1=b(1:n,:);
            b2=b(n+1:2*n,:);
            b3=b(2*n+1:3*n,:);
            b4=b(3*n+1:end,:);
            btake1=zeros(n,lp);
            btake2=btake1;
            btake3=btake1;
            btake4=btake1;
            for nn=1:n
                movenn=toAge(nn);
                %
                btakenn=zeros(1,lp);
                bnn=b1(nn,:);
                findnn=find(bnn>0);
                popsnn=bnn(findnn);
                fromnn=randsample(1:NNbar(nn),movenn);
                [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
                btakenn(findnn)=histVec;
                btake1(nn,:)=btakenn;
                %
                btakenn=zeros(1,lp);
                bnn=b2(nn,:);
                findnn=find(bnn>0);
                popsnn=bnn(findnn);
                fromnn=randsample(1:NNbar(n+nn),movenn);
                [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
                btakenn(findnn)=histVec;
                btake2(nn,:)=btakenn;
                %
                btakenn=zeros(1,lp);
                bnn=b3(nn,:);
                findnn=find(bnn>0);
                popsnn=bnn(findnn);
                fromnn=randsample(1:NNbar(2*n+nn),movenn);
                [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
                btakenn(findnn)=histVec;
                btake3(nn,:)=btakenn;
                %
                btakenn=zeros(1,lp);
                bnn=b4(nn,:);
                findnn=find(bnn>0);
                popsnn=bnn(findnn);
                fromnn=randsample(1:NNbar(3*n+nn),movenn);
                [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
                btakenn(findnn)=histVec;
                btake4(nn,:)=btakenn;
            end
        else
            btake1=b(1:n,:).*repmat(vover1,1,lp);%%Diff for age
            btake2=b(n+1:2*n,:).*repmat(vover2,1,lp);
            btake3=b(2*n+1:3*n,:).*repmat(vover3,1,lp);
            btake4=b(3*n+1:end,:).*repmat(vover4,1,lp);
        end

        b=b-[btake1;btake2;btake3;btake4]+[zeros(n,lp);btake1;btake2;btake3];
        b(1:n,1)=b(1:n,1)+sum(btake4,2);
        if min(min(min(b)))<0
            error('Shit the bed')
        end
        %%
        %New loop (for age structure):
        if t<tauend
            if solvetype==3
                error('Code not written yet')
            else
                %Define NNbar, NN etc
                %[gamma,NN,n,ntot,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Q,0,0);
                %Don't change D's - cheating?; total pops in each location not
                %changing
                NNall=circshift(NNall,1,2);
                Na1=sum(NNall(:,1:5),2); Na2=sum(NNall(:,6:20),2); Na3=sum(NNall(:,21:65),2); Na4=sum(NNall(:,66:end),2);
                NNbar=[Na1,Na2,Na3,Na4];
                vover1=NNall(:,5)./Na1; vover2=NNall(:,20)./Na2; vover3=NNall(:,65)./Na3; vover4=NNall(:,end)./Na4; 
                vover1(Na1==0)=1; vover2(Na2==0)=1; vover3(Na3==0)=1; vover4(Na4==0)=1;
                NNbar=reshape(NNbar,ntot,1);
            end
        end
    end
    S0=b(:,1);
    Z0=sum(b(:,2:end),2);%*(1-mu);
end
end
%%
%}
f=A1;
g=A2;
end
%%
function f=solveZi(Zi,Z0,beta,gamma,D,Nages,addbit)
lZi=log(Nages-Zi); lZi(Zi>1)=NaN;
f=lZi-log(Nages-Z0)+beta/gamma*D*(Zi-Z0)+addbit;
end