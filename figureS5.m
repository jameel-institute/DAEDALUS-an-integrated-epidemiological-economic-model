close all;
clear all;

addpath('../');
addpath('../Output Data/GDP-63/');

load('UK63.mat','data');
load('B2.mat');

numInt=length(xoptim)/length(data.G);
tvec=[-61.3373  87.6062  245 306  367  426];
t=[1:1:tvec(end)];

sdfact=0.05;
numIter=1000;
Hf=zeros(tvec(end),numIter);
Hc=zeros(tvec(end),numIter);
He=zeros(tvec(end),numIter);
rng default;%for reproducibility

%%

for i=1:numIter
    
    ddataxi=data;
    
    ddataxi.comm=normrnd(ddataxi.comm,sdfact*ddataxi.comm);
    
    ddataxi.schoolA1=normrnd(ddataxi.schoolA1,sdfact*ddataxi.schoolA1);
    ddataxi.schoolA2=normrnd(ddataxi.schoolA2,sdfact*ddataxi.schoolA2);
    
    ddataxi.travelA3=normrnd(ddataxi.travelA3,sdfact*ddataxi.travelA3);
    ddataxi.hospA2=normrnd(ddataxi.hospA2,sdfact*ddataxi.hospA2);
    ddataxi.hospA3=normrnd(ddataxi.hospA3,sdfact*ddataxi.hospA3);
    ddataxi.hospA4=normrnd(ddataxi.hospA4,sdfact*ddataxi.hospA4);
    
    ddataxi.B=normrnd(ddataxi.B,sdfact*ddataxi.B);
    ddataxi.C=normrnd(ddataxi.C,sdfact*ddataxi.C);
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(ddataxi,numInt);%,inp);
    pr.sw=0;%switching off
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,ddataxi);
    
    tt=f(:,1);
    hh=f(:,4);
    Hf(:,i)=interp1(tt,hh,t);
    
end

Hf=prctile(Hf,[5,50,95],2);
f1=Hf(:,1);
f2=Hf(:,2);
f3=Hf(:,3);

%%

for i=1:numIter
    
    ddataxi=data;
    
    ddataxi.comm=normrnd(ddataxi.comm,sdfact*ddataxi.comm);
    
    ddataxi.schoolA1=normrnd(ddataxi.schoolA1,sdfact*ddataxi.schoolA1);
    ddataxi.schoolA2=normrnd(ddataxi.schoolA2,sdfact*ddataxi.schoolA2);
    
    ddataxi.travelA3=normrnd(ddataxi.travelA3,sdfact*ddataxi.travelA3);
    ddataxi.hospA2=normrnd(ddataxi.hospA2,sdfact*ddataxi.hospA2);
    ddataxi.hospA3=normrnd(ddataxi.hospA3,sdfact*ddataxi.hospA3);
    ddataxi.hospA4=normrnd(ddataxi.hospA4,sdfact*ddataxi.hospA4);
    
    %ddataxi.B=normrnd(ddataxi.B,sdfact*ddataxi.B);
    %ddataxi.C=normrnd(ddataxi.C,sdfact*ddataxi.C);
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(ddataxi,numInt);%,inp);
    pr.sw=0;%switching off
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,ddataxi);
    
    tt=f(:,1);
    hh=f(:,4);
    Hc(:,i)=interp1(tt,hh,t);
    
end

Hc=prctile(Hc,[5,50,95],2);
c1=Hc(:,1);
c2=Hc(:,2);
c3=Hc(:,3);

%%

for i=1:numIter
    
    ddataxi=data;
    
    %ddataxi.comm=normrnd(ddataxi.comm,sdfact*ddataxi.comm);
    
    ddataxi.schoolA1=normrnd(ddataxi.schoolA1,sdfact*ddataxi.schoolA1);
    ddataxi.schoolA2=normrnd(ddataxi.schoolA2,sdfact*ddataxi.schoolA2);
    
    %ddataxi.travelA3=normrnd(ddataxi.travelA3,sdfact*ddataxi.travelA3);
    %ddataxi.hospA2=normrnd(ddataxi.hospA2,sdfact*ddataxi.hospA2);
    %ddataxi.hospA3=normrnd(ddataxi.hospA3,sdfact*ddataxi.hospA3);
    %ddataxi.hospA4=normrnd(ddataxi.hospA4,sdfact*ddataxi.hospA4);
    
    %ddataxi.B=normrnd(ddataxi.B,sdfact*ddataxi.B);
    %ddataxi.C=normrnd(ddataxi.C,sdfact*ddataxi.C);
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(ddataxi,numInt);%,inp);
    pr.sw=0;%switching off
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,ddataxi);
    
    tt=f(:,1);
    hh=f(:,4);
    He(:,i)=interp1(tt,hh,t);
    
end

He=prctile(He,[5,50,95],2);
e1=He(:,1);
e2=He(:,2);
e3=He(:,3);

%%

numPeriods=length(tvec)-1;
numSectors=length(data.G);
lt=length(tvec);
dodiff=1;

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
hold on;
lw=2.0;
thresh=[12000,18000,24000];
numThresh=length(thresh);

maxY=48000;
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=2
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end

t2=[t,fliplr(t)];
inBetween=[f1',fliplr(f3')];
fill(t2,inBetween,'red','facealpha',.2);
plot(t,f1','-','linewidth',lw/2,'color','red');
plot(t,f3','-','linewidth',lw/2,'color','red');
h1=plot(t,f2','-','linewidth',lw,'color','red');

inBetween=[c1',fliplr(c3')];
fill(t2,inBetween,'black','facealpha',.2);
plot(t,c1','-','linewidth',lw/2,'color','black');
plot(t,c3','-','linewidth',lw/2,'color','black');
h2=plot(t,c2','--','linewidth',lw,'color','black');

points=tvec+10;
pointsy=.93*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',15);
text(tvec(2)+5,pointsy,'LD','fontsize',15);
for i=3:lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',15)
end
xlim([0,tvec(end)]);
ylim([0,maxY]);
axis square;
xlabel('Time','FontSize',fs);
ylabel('Hospital Occupancy','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-20 0 0]);
%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[1,32,61,92,122,153,183,214,245,275,306,336,367,398]);
set(gca,'ytick',[0:6000:maxY]);
set(gca,'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'});
if numPeriods==5
	xlabels2=({'PRE','LD','1','2','3'});
elseif numPeriods==8
    %xlabels2=({'Jan','Mar 26th','Sep','Nov','Jan'});
    xlabels2=({'PRE','LD','1','2','3','4','5','6'});
else
    error('Data missing for nunmPeriods')
end
xtickangle(45);
ax = gca;
ax.YAxis.Exponent = 3;
box on;
grid on;
grid minor;
%legend([hh1,hh2],'Inc.','Hosp. occ.','location','west')
%legend([hh1,hh2,hh3,hh4,hh5],'Incidence','I','H','D','V','location','northwest');
%legend([hh0,hh4,hh1,hh2,hh3],'LDA','FO','A (12,000)','A (18,000)','A (24,000)','Position',[-0.245 0.24 1 1]);
legend([h1,h2],'Full','Community','Position',[-0.239 0.288 1 1]);

% pointsx=385;
% pointsy=10^3*[2.5,13.5,19.5,25.5,82];
% txt=[hlda,ha1,ha2,ha3,hfo];
% for i=1:5
%     text(pointsx,pointsy(i),['£' num2str(txt(i)) 'bn'],'fontsize',12);
% end

hold off;

%%

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
hold on;
lw=2.0;
thresh=[12000,18000,24000];
numThresh=length(thresh);

maxY=48000;
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=2
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end

t2=[t,fliplr(t)];
inBetween=[f1',fliplr(f3')];
fill(t2,inBetween,'red','facealpha',.2);
plot(t,f1','-','linewidth',lw/2,'color','red');
plot(t,f3','-','linewidth',lw/2,'color','red');
h1=plot(t,f2','-','linewidth',lw,'color','red');

inBetween=[e1',fliplr(e3')];
fill(t2,inBetween,'black','facealpha',.2);
plot(t,e1','-','linewidth',lw/2,'color','black');
plot(t,e3','-','linewidth',lw/2,'color','black');
h2=plot(t,e2','--','linewidth',lw,'color','black');

points=tvec+10;
pointsy=.93*maxY;
txt={'1','2','3','4','5','6'};
text(10,pointsy,'PRE','fontsize',15);
text(tvec(2)+5,pointsy,'LD','fontsize',15);
for i=3:lt-1
    text(points(i),pointsy,txt{i-2},'fontsize',15)
end
xlim([0,tvec(end)]);
ylim([0,maxY]);
axis square;
xlabel('Time','FontSize',fs);
ylabel('Hospital Occupancy','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-20 0 0]);
%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[1,32,61,92,122,153,183,214,245,275,306,336,367,398]);
set(gca,'ytick',[0:6000:maxY]);
set(gca,'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'});
if numPeriods==5
	xlabels2=({'PRE','LD','1','2','3'});
elseif numPeriods==8
    %xlabels2=({'Jan','Mar 26th','Sep','Nov','Jan'});
    xlabels2=({'PRE','LD','1','2','3','4','5','6'});
else
    error('Data missing for nunmPeriods')
end
xtickangle(45);
ax = gca;
ax.YAxis.Exponent = 3;
box on;
grid on;
grid minor;
%legend([hh1,hh2],'Inc.','Hosp. occ.','location','west')
%legend([hh1,hh2,hh3,hh4,hh5],'Incidence','I','H','D','V','location','northwest');
%legend([hh0,hh4,hh1,hh2,hh3],'LDA','FO','A (12,000)','A (18,000)','A (24,000)','Position',[-0.245 0.24 1 1]);
legend([h1,h2],'Full','Education','Position',[-0.247 0.288 1 1]);

% pointsx=385;
% pointsy=10^3*[2.5,13.5,19.5,25.5,82];
% txt=[hlda,ha1,ha2,ha3,hfo];
% for i=1:5
%     text(pointsx,pointsy(i),['£' num2str(txt(i)) 'bn'],'fontsize',12);
% end

hold off;