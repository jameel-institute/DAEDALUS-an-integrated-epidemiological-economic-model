function [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt)%,R0,del)%,inp)

%data.NNs - column vector of population
%Possible generalsiation to within-sector heterogeneity - one column per subsector

%% POPULATION PARAMETERS:

%Population Density
[n,na]=size(data.NNs);
ntot=n*na;
NN=sum(data.NNs,2);
NNbar=reshape(data.NNs,ntot,1);
NNrep=repmat(NN,na,1);

%Urban and Rural?
%urbrur=0;%Turn in to @home vs @work7/291
%{
if urbrur==1
    hvec=kron(hvec,ones(n,1));
    muvec=kron(muvec,ones(n,1));
    %Kkron=[1,.05;.05,.75];%Sonoma
    Kkron=[1,.05;.05,.85];%Mendocino
else
    Kkron=1;
end
%}
%{
C=[.3827    0.9115    0.0419;
    1.2062    7.3538    0.1235;
    0.1459    0.4810    0.1860];
%}
%C=eye(na);

%Age-Sector Breakdown
lc=4;
adInd=3;
lx=length(data.NNs)-lc;
D=heMakeDs(NN,ones(lx,1),data,0);
Dout=D;

%K=heMakeDs(NN,eye(10));
%K=rand(n);
%K=normr(K);
%D=kron(C,K);

pr=struct;%from Knock et al. (2020) unless indicated otherwise
pr.a=1;%0.59;

%% DISEASE PARAMETERS:

%Transmission
pr.R0=2.7500;%R0;%from fitting
pr.red=0.58;%from Byambasuren et al. (2020)

%Latency and Onset
Text=4.6;
pr.sigma=1/Text;
%Tonset=1;
%pr.omega=1/Tonset;

%Case Pathways
pr.p1=0.6;
[ph,pd,~,~]=heParamsAge(data,pr.p1);
% pili=[0.4670,0.4786,0.6590,0.7252]';
% pili=0.505*pili;
% pili=[repmat(pili(adInd),lx,1);pili];
% pr.p2=pili;
% ph=[0.0500    0.0022    0.0350    0.5235]';
ph=[repmat(ph(adInd),lx,1);ph];
% pd=[0.0103    0.0078    0.0361    0.1555]';
pd=[repmat(pd(adInd),lx,1);pd];
% pdeath=0.39;
% pd=pd/sum(pd);
% pd=pdeath*pd;%new data: 39% of hospital cases result in death

%Hospitalisation and Death Rates (proportion and rate combined)
Tsh=4.5;%from Global
pr.h=ph/Tsh;
Thd=11.1;%from Pablo
pr.mu=pd/Thd;

%Recovery Rates (proportion and rate combined)
Ta=2.1;%asymptomatic
pr.g1=1/Ta;
% pr.gX=1/2.1;%mild
Ts=4;%symptomatic
pr.g2=(1-ph)/Ts;
Threc=11.1;%hospitalised%from Pablo
pr.g3=(1-pd)/Threc;

%Immunity Loss
Ti=150;%from RTM
pr.nu=0;%1/Ti;

%%

Deff=Dout.*repmat(NNbar,1,n*na)./repmat(NNbar',n*na,1);
onesn=ones(ntot,1);

F=zeros(3*ntot,3*ntot);
%F=zeros(4*ntot,4*ntot);
F(1:ntot,ntot+1:end)=[pr.red*Deff,Deff];
%F(1:ntot,ntot+1:end)=[pr.red*Deff,repmat(Deff,1,2)];

vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       (pr.g2+pr.h).*onesn];%g2 and h are vectors
%vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       pr.g2.*onesn;       (pr.gX+pr.h).*onesn];%gX and h are vectors
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=    diag(-(1-pr.p1) .*pr.sigma  .*onesn);
V(2*ntot+1:3*ntot,1:ntot)=  diag(-pr.p1     .*pr.sigma  .*onesn);
%V(2*ntot+1:3*ntot,1:ntot)=  diag(-pr.p1     .*(1-pr.p2) .*pr.sigma  .*onesn);
%V(3*ntot+1:4*ntot,1:ntot)=  diag(-pr.p1     .*pr.p2     .*pr.sigma  .*onesn);

GD=F/V;

%Ceff=kron(C,Ckron);%Urb/rural mixing here
%F(1:ntot,1:ntot)=Deff;
%vvec=kron([pr.sigma;pr.g1;pr.g2;pr.gX],ones(ntot,1));
%vvec(end-ntot+1:end)=vvec(end-ntot+1:end)+pr.h;

%{
%HE:
ntot=n*na;
F=zeros(7*ntot,7*ntot);
F(1:ntot,1:ntot)=Deff;%ntot+1:end)=[repmat(2/3*Deff,1,2),repmat(Deff,1,4)];
onesn=ones(ntot,1);
vvec=[pr.sigma*onesn;pr.g1*onesn;pr.omega*onesn;pr.g2*onesn;pr.g2*onesn;pr.h*onesn;pr.h*onesn];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*onesn);
V(3*ntot+1:4*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p3).*(1-pr.p2));
V(4*ntot+1:5*ntot,2*ntot+1:3*ntot)=diag(-pr.p3.*(1-pr.p2));
V(5*ntot+1:6*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p4).*pr.p2);
V(6*ntot+1:7*ntot,2*ntot+1:3*ntot)=diag(-pr.p4.*pr.p2);
GD=F/V;
%}

d=eigs(GD,1);%largest in magnitude (+/-) 
R0a=max(d); 
beta=pr.R0/R0a;%beta scales actual R0 to fitted R0

%% PREPAREDNESS PARAMETERS:

%Mitigation Time
pr.Tm=72.240625860903876;%121;%sector closures (x), home-working (wfhAv), NPIs (betamod) and self-isolation (p3,p4,p5) implemented

%Adherence to NPIs
% pmod=0.5408;%lockdown delta%from fitting
% pmodIn=0.6;%post-lockdown delta
pmod=0.5650;%0.6227;%lockdown delta%from fitting
pmodIn=0.72;%0.79;%post-lockdown delta
pr.betamod=[1,pmod,linspace(pmodIn,1,numInt)];%[1,del];

%Self-Isolation and Quarantine
pr.p3=0.00;%proportion of asymptomatic self-isolating
pr.p4=0.00;%proportion of symptomatic self-isolating
%pr.p5=0.00;%proportion of severe self-isolating
%
pr.odds=0;%asymptomatic quarantining rate
pr.q1=0;%mild quarantining rate
pr.q2=0;%severe quarantining rate
% pr.g4=1/(1/pr.g2-1/pr.q1);%asymptomatic/mild quarantining recovery rate
% pr.qnew=0;
% if pr.q2>0
%     pr.g4X=1/(1./pr.gX+1./pr.q2);%Vector
% else
%     pr.g4X=pr.q2*0;
% end

%Hospital Capacity
pr.Hmax=18000;
pr.mu_oc=pr.mu;%1.19*pr.mu;%Wilde et al. (2021)
pr.g3_oc=pr.g3;

%% VACCINATION PARAMETERS:

vx=struct;%from RTM unless indicated otherwise

%Vaccine 1
vx.hrv1=    1/21;   %time to develop v-acquired immunity (AstraZeneca)
vx.scv1=    0;      %infection-blocking efficacy
vx.p1v1=    0;      %disease-blocking efficacy          
% vx.hv1=   0;      %severe-disease-blocking efficacy% vx.p2v1=0;
vx.trv1=    0;      %transmission-blocking efficacy
vx.nuv1=    1/365;  %duration of v-acquired immunity 

% %Vaccine 2
% vx.hrv2=    1/21;   %time to develop v-acquired immunity (Pfizer)
% vx.scv2=    0;      %infection-blocking efficacy
% vx.p1v2=    0;      %disease-blocking efficacy          
% % vx.hv2=   0;      %severe-disease-blocking efficacy% vx.p2v2=0;
% vx.trv2=    0;      %transmission-blocking efficacy
% vx.nuv2=    1/365;  %duration of v-acquired immunity 

%Population
Npop=   data.Npop;
NNage=  [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
prop50= sum(Npop(11:13))/sum(Npop(5:13));%proportion of working-age adults over 50
%propCare=398840/vx.NNage(4);

%Uptake
upUnder50=  0.75;
upOver50=   1.00;
%upCare=    0.95;

%Rollout
tpoints=[336,367,457];%rollout periods: 1st December, 1st January, 1st April
tint=   diff(tpoints);%duration of rollout periods

%Period 1 - 2 million doses to over 65s
vx.startp1= tpoints(1);
nump1=      [0;0;0;0];%[0;0;0;2*10^6];%by age group only
vx.aratep1=  nump1/tint(1);

%Period 2 - doses for all over 50s (based on 100% uptake)
vx.startp2= tpoints(2);
nump2=      [0;0;0;0];%[0;0;prop50*NNage(3);NNage(4)-2*10^6];%by age group only
vx.aratep2=  nump2/tint(2);%to be split across all economic sectors in heSimCovid19vax.m

% Calculate which group exhausts first (accounting for uptake):
% toVax=NNnext.*[0;0;.75;.95*propCare+.75*(1-propCare)];%Accounting for uptake
% T=toVax./vx.rate2age;%Pop to vaccinate/vax rate %Assumes NNage(4)>2e6
% T(vx.rate2age==0)=0;
% [tnext,which]=min(T);
% T=[0;0;68.25;68.8081];

%Period 3 - doses for all remaining working-age
vx.startp3= tpoints(2)+69;%start time for Period 3 based on block above (not from tpoints as that assumes 100% uptake)
ratep2=     sum(vx.aratep2);%rate for remainder of the rollout
vx.aratep3=  [0;0;ratep2;0];%to be split across all economic sectors in heSimCovid19vax.m

%End of Rollout
toVax=upUnder50*(1-prop50)*NNage(3);%accounting for uptake
T=toVax/ratep2;%remaining population to vaccinate/vaccination rate
vx.end=vx.startp3+T;

end