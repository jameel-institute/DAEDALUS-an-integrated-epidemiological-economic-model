function poptim=fitEpi20(ydata,xminData)
%%

load('UK63.mat','data');

xmin=table2array(xminData(:,5:24));
xmin=min(xmin/100,1);%March values mixed up

% xmins(:,4)= (113001.035*xmins(:,4)    +37981.16*xmins(:,5))    /(113001.035    +37981.16);%IO table, UK, 2016
% xmins(:,13)=(249098.011*xmins(:,13)   +164492.031*xmins(:,14)) /(249098.011    +164492.031);%IO table, UK, 2016
% xmins(:,18)=(47439.008*xmins(:,18)    +40305.011*xmins(:,19))  /(47439.008     +40305.011);%IO table, UK, 2016
% xmins(:,19)=[];
% xmins(:,14)=[];
% xmins(:,5)=[];
% X=zeros(11,36);
% repCol=[1,3,16,1,1,1,1,1,3,1,1,1,1,1,1,1,1];
X=zeros(12,63);
repCol=[3,1,19,1,2,1,3,5,1,4,3,1,5,4,1,1,2,2,3,1];

repColSum=[0,cumsum(repCol)];
for i=1:length(repCol)
    X(:,repColSum(i)+1:repColSum(i+1))=repmat(xmin(:,i),1,repCol(i));
end

X=X(6:9,:);%May-August incl.
%X=X(1:7,:);%May-August incl.
%X=repmat(data.xmin,4,1);%model approximation

%
X=reshape(X',63*4,1);
%X=reshape(X',36*7,1);

ydata=table2array(ydata);%England only, 20th March-4th February incl.
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
xdata=(80:244);%20th March-31st August incl.%days of year
ydata=ydata(1:size(xdata,2));

ub=[0,      100,    2.75,   1.00*ones(1,2)];
x0=[-50,    90,     2.75,   0.50*ones(1,2)];
lb=[-100,   60,     1.01,   0.00*ones(1,2)];%t0, tld, R0, deltas

%%

fun=@(params,xdata)sim2fit(params,data,xdata,X);

rng default;%for reproducibility
options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1000000);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata,'lb',lb,'ub',ub,'options',options);
ms=MultiStart;
poptim=run(ms,problem,1);

% [poptim,~,R,~,~,~,J]=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
% CovB=(1/(J'*J))*((R'*R)./(numel(ydata) - numel(poptim)));
% [poptim,R,J,CovB,~,~]=nlinfit(xdata,ydata,fun,x0);
% conf=nlparci(poptim,R,'jacobian',J);
% [Ypred,delta]=nlpredci(fun,xdata,poptim,R,'Covar',CovB);
% f=poptim;
% g=conf;

ymod=sim2fit(poptim,data,[round(poptim(1)):xdata(end)],X);

f=figure('Units','centimeters','Position',[0 0 20 20]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',15);
hold on;
bar(xdata,ydata);
plot([round(poptim(1)):xdata(end)],ymod,'linewidth',2.5,'color','red');
xlim([xdata(1),xdata(end)]);
axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Time');
ylabel('Hospital Occupancy');
title('Model Fit');

end


function f=sim2fit(params,data,xdata,Xfit)
    
    numInt=length(Xfit)/length(data.G);
    
    tvec=[params(1),params(2),122,153,183,214,245];
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,params(3),[params(4)*ones(1,2),params(5)*ones(1,3)]);
    pr.sw=0;
    Wfit=Xfit.^(1/pr.a);
    
    [simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec,0,data);

    % if simu(1,1)>1
    %     simu=[zeros(simu(1,1)-1,1);simu(:,4)];%???
    % else
    %     simu=simu;
    % end

    t=simu(:,1)';
    h=simu(:,4)';

    f=interp1(t,h,xdata); 

end