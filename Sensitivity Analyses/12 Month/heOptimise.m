function [xoptim,yoptim,exitflag]=heOptimise

%Hosptial capacity defined in hePrepCovid19.m

%% INITIALISATION:

load('UK63.mat','data');

numSect=length(data.G);
numInt=6;
lx=numInt*numSect;
monthPeriod=2;%intervention intervals (1 for 6x1, 2 for 3x2)

[pr,vx,NN,n,ntot,na,NNbar,NNrep,Din,beta]=hePrepCovid19(data,numInt);
pr.sw=0;%switching off

% tvec=[-61.3373   87.6062  245 306  367  426];
% tvec=[-60.7399   86.8162  245 306  367  426];%ChLessSus fit
% tvec=[-61.3373   87.6062  245 275  306  336  367  398  426];%6x1
tvec=[-61.3373   87.6062  245 306  367  426  487  548  610];%6x2
% tvec=[1  linspace(pr.Tm,  365*2+1,  numInt+1)];

scenB=1;%scenario A or B

X0in=0;%option to input X0%if numPeriod=1, code assumes that X0 is a 3x2 solution
if X0in==1
    load('IC.mat','xoptim');
end

fileName='xoptim.mat';%output file name

%% CONSTRAINTS:

%Epidemiological constraints (H_max, R_end)
hospThresh=[pr.Hmax,1];

%Economic constraints (over whole period)
G=data.G*monthPeriod;
Z2=repmat(G,1,numInt);
%b2=zeros(numSect,1);
xmn=min(data.xmin',1);
xmn=0.8*xmn;
if numSect==36
    xmn(34)=data.xmin(34);
elseif numSect==63
    xmn(56)=data.xmin(56);
else
    error('Unknown economic configuration!');    
end
b2=min(numInt*G*xmn,0);%multiplied by monthPeriod above

%Z2=[repmat(G,1,numInt);repmat(-G,1,numInt)];
%b2=[b+1e-6;zeros(numSect,1)];
%b2=[b;zeros(numSect,1)]+4000;

xmin=repmat(data.xmin',numInt,1);%lockdown A

ub=ones(lx,1);
ub(xmin>1)=xmin(xmin>1);

xmin(xmin>1)=1;

lb=0.8*xmin;
if numSect==36
    lb(34:numSect:end)=xmin(34);
    if scenB==1
        lb(33:numSect:end)=0.8;
    end
elseif numSect==63
    lb(56:numSect:end)=xmin(56);
    if scenB==1
        lb(55:numSect:end)=0.8;
    end
else
    error('Unknown economic configuration!');
end

if X0in==0
    X0=xmin;
    if scenB==1
        if numSect==36
            X0(33:numSect:end)=0.8;
        elseif numSect==63
            X0(55:numSect:end)=0.8;
        else
            error('Unknown economic configuration!');
        end
    end
    %options=optimoptions('linprog','OptimalityTolerance',1e-10);
    X0=linprog(ones(1,size(Z2,2)),Z2,b2,[],[],X0,X0+0.90*(ub-X0));
else
    X0=xoptim;
    %     if monthPeriod==1
    %         X0=[X0(1:63);X0(1:63);X0(64:126);X0(64:126);X0(127:189);X0(127:189)];
    %     end
end

objFun=data.obj*monthPeriod;
objFun=repmat(objFun,numInt,1);
% objFunemp=NNsector(1:numSect);
% objFunemp=repmat(objFunemp,numInt,1);

fun1=@(Xit)econGDP(objFun,Xit);
%fun1=@(Xit)employment(objFunemp,Xit);
nonlcon=@(Xit)epiConstraint(pr,vx,n,ntot,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,data);

%%

%{
options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'algorithm','interior-point');%'sqp' %'interior-point'
[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%[xoptim,yoptim,exitflag]=fmincon(fun1,X0,Z2,b2,[],[],lb,ub,[],options);
%}
%
%{
rng default % For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
ms=MultiStart('StartPointsToRun','bounds');%,'Display','iter','MaxTime',3600);
%rs=RandomStartPointSet('NumStartPoints',3);
%allpts={rs};
%[x fval eflag output manymins]=run(ms,problem,allpts);
%
%numPoints=12;
%sp=repmat(lb',numPoints,1)+repmat(ub'-lb',numPoints,1).*rand(numPoints,lx);
%sp=CustomStartPointSet(sp);
[xoptim,yoptim,exitflag,output,manymins]=run(ms,problem);%,sp);
%exitflag='Not relevant for ms';
%}
%
%{
rng default
options=optimoptions('patternsearch','UseParallel',true,'UseVectorized', false,'UseCompletePoll',true);%,'InitialPoints',ips);%'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);
%,'ConstraintTolerance',1);
[xoptim,yoptim,exitflag]=patternsearch(fun1,X0,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
%
%{
rng default
%ips=struct;
%ips.X0=[.9*X0';.8*X0'+.5*(ub'-X0')];
np=23;
xmin=repmat(.9*X0',np,1);
xmax=repmat(.1*X0'+.5*(ub'-X0'),np,1);
ips=xmin+rand(np,378/monthPeriod).*xmax;
ips=[X0';ips];
fun2=@(Xit)econGDPga(objFun,Xit);
optimoptions('patternsearch','UseVectorized',true);%,'UseCompletePoll',true);
options=optimoptions('paretosearch','InitialPoints',ips);%'UseParallel',false,'UseVectorized', true,'InitialPoints',ips);%,'ParetoSetSize',1e4 ,'UseCompletePoll',true
[xoptim,yoptim,exitflag]=paretosearch(fun2,378/monthPeriod,Z2,b2,[],[],lb,ub,nonlcon,options);
%}
%
%{
fun2=@(Xit)econGDPga(objFun,Xit);
rng default
options=optimoptions('ga','UseParallel',true,'UseVectorized',false,'InitialPopulationRange',[X0';X0'+.1*(ub'-X0')]);
%opts.InitialPopulationRange=[lb'; X0'];%2 rows, nvar columns - lb;ub
[xoptim,yoptim,exitflag]=ga(fun2,378/monthPeriod,Z2,b2,[],[],lb,ub,nonlcon,options);
xoptim=xoptim';
%}

%% OPTIMISATION:

rng default;%for reproducibility
options=optimoptions(@fmincon,'UseParallel',true,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);%,'algorithm','interior-point');
problem=createOptimProblem('fmincon','x0',X0,'objective',fun1,'Aineq',Z2,'bineq',b2,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
gs=GlobalSearch;
[xoptim,yoptim,exitflag]=run(gs,problem);
save(fileName,'xoptim','yoptim','exitflag');

end

%% FUNCTIONS:

function f=econGDP(obj,Xit)

    f=-sum(obj.*Xit);

end

function f=employment(obj,Xit)

    f=-sum(obj.*Xit);

end

function [c,cex]=epiConstraint(pr,vx,n,ntot,na,NN,NNbar,NNrep,Din,beta,Xit,tvec,hospThresh,data)
    
    Wit=Xit.^(1/pr.a);
    [~,h,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Din,beta,Wit,tvec,0,data);
    c=h-hospThresh;
    cex=[];

end

%%

%{
%Econ constrained by month:
ball=repmat(b',numInt,1);
kkron=eye(numInt-1);
kkron(end+1,end+1)=1/lastInt;%Account for last interval being longer
Gall=sparse(kron(kkron,G));
lx=size(Gall,1);
%Gvec=diag(Gall);
Z2=[Gall;-Gall];%<0,b
b2=[ball/numInt;zeros(lx,1)];
%
lb=zeros(lx,1);
ub=ones(lx,1);
X0=.1*rand(lx,1);
%}
%{
%G(:,44)=G(:,44)+G(:,45);
%G(44,:)=G(44,:)+G(45,:);
G(45,:)=[];
G(:,45)=[];
%b(44)=b(44)+b(45);
b(45)=[];
%objFun(44)=objFun(44)+objFun(45);
objFun(45)=[];
%
%data64.xmin=zeros(1,63);
b(2)=0;%b(1);
%
objFun(2)=objFun(1);
g1r=G(1,:);
gc1=G(:,1);
G(2,:)=g1r;
G(:,2)=gc1;
%}
%
%X0=repmat(data.xmin',numInt,1);%zeros(lx,1);
%X0=[ones(2*63,1);repmat(data.xmin',4,1)];
%xunder=data.xmin'+.4*(1-data.xmin');
%X0=[repmat(xunder,2,1);repmat(data.xmin',4,1)];
%X0(end-62:end)=data.xmin';
%
%{
%Schools open fully/to given degree:
xlb=repmat(data.xmin',numInt,1);
xlb((0:5)*63+55,1)=1;
xlb(xlb>1)=1;
lb=xlb;
X0=xlb;
%}
%X0(X0<.7)=.7;
%X0(end-62:end)=data.xmin';
%
%X0(2*63+55:63:end)=.8;
%Start-point based on existing output:
%X0=xoptim;
%xmin64=repmat(data.xmin',6,1);
%X0=xmin64+.8*(xoptim-xmin64);
%X0(5*63+1:end)=xmin64(5*63+1:end);
%X0=xmin64+kron([.5;.4;.3;.2;.1;0],ones(63,1)).*(xoptim-xmin64);
%Schox0(x0<.7)=.7;ols open:
%X0(4*63+55,1)=1; X0(5*63+55,1)=1;