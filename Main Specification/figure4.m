close all;
clear all;

addpath('../');
load('UK63.mat','data');
addpath('../Output Data/GDP-63/');
load('A1.mat');
numInt=length(xoptim)/length(data.G);

[pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt);%,inp);
pr.sw=0;%switching off
tvec=[-61.3373   87.6062  245 306  367  426];

xoptim=repmat(data.xmin',3,1);
[flda,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
hlda=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

load('A1.mat')
[fa1,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
ha1=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

load('A2.mat');
[fa2,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
ha2=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

load('A3.mat');
[fa3,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
ha3=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

xoptim=ones(189,1);
[ffo,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
hfo=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

%%

numPeriods=length(tvec)-1;
numSectors=length(data.G);
lt=length(tvec);
dodiff=1;

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
lw=2.5;
cmap=lines(2);
thresh=[12000,18000,24000];
numThresh=length(thresh);

tlda=flda(:,1); 
hoslda=flda(:,4);

ta1=fa1(:,1);
hosa1=fa1(:,4);

ta2=fa2(:,1);
hosa2=fa2(:,4);

ta3=fa3(:,1);
hosa3=fa3(:,4);

tfo=ffo(:,1);
hosfo=ffo(:,4);

maxY=96000;
%
hold on;
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=1:numThresh
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end

hh0=plot(tlda,hoslda,'-','linewidth',lw,'color',[0.5 0.5 0.5]);
hh4=plot(tfo,hosfo,':','linewidth',lw,'color',[0.5 0.5 0.5]);
hh1=plot(ta1,hosa1,'-','linewidth',lw,'color','blue');
hh2=plot(ta2,hosa2,'-','linewidth',lw,'color','red');
hh3=plot(ta3,hosa3,'-','linewidth',lw,'color',[0.9290, 0.6940, 0.1250]);
hh5=plot(tfo(1:250),hosfo(1:250),'-','linewidth',lw,'color','black');

% points=tvec+10;
% pointsy=.93*maxY;
% txt={'1','2','3','4','5','6'};
% text(10,pointsy,'PRE','fontsize',15);
% text(tvec(2)+5,pointsy,'LD','fontsize',15);
% for i=3:lt-1
%     text(points(i),pointsy,txt{i-2},'fontsize',15)
% end
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
set(gca,'ytick',[0:12000:96000]);
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
legend([hh0,hh4,hh1,hh2,hh3],'LDA','FO','A (12,000)','A (18,000)','A (24,000)','Position',[-0.245 0.325 1 1]);

pointsx=385;
pointsy=10^3*[2.5,13.5,19.5,25.5,82];
txt=[hlda,ha1,ha2,ha3,hfo];
for i=1:5
    text(pointsx,pointsy(i),['£' num2str(txt(i)) 'bn'],'fontsize',12);
end

hold off;

%%

xoptim=repmat(data.xmin',3,1);
xoptim(55:63:end)=0.80;
[flda,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
hlda=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

load('B1.mat')
[fa1,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
ha1=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

load('B2.mat');
[fa2,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
ha2=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

load('B3.mat');
[fa3,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
ha3=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

xoptim=ones(189,1);
[ffo,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
hfo=round(sum(xoptim.*repmat((6/numInt)*data.obj,numInt,1))/1000);

%%

numPeriods=length(tvec)-1;
numSectors=length(data.G);
lt=length(tvec);
dodiff=1;

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
lw=2.5;
cmap=lines(2);
thresh=[12000,18000,24000];
numThresh=length(thresh);

tlda=flda(:,1); 
hoslda=flda(:,4);

ta1=fa1(:,1);
hosa1=fa1(:,4);

ta2=fa2(:,1);
hosa2=fa2(:,4);

ta3=fa3(:,1);
hosa3=fa3(:,4);

tfo=ffo(:,1);
hosfo=ffo(:,4);

maxY=96000;
%
hold on;
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=1:numThresh
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end

hh0=plot(tlda,hoslda,'-','linewidth',lw,'color',[0.5 0.5 0.5]);
hh4=plot(tfo,hosfo,':','linewidth',lw,'color',[0.5 0.5 0.5]);
hh1=plot(ta1,hosa1,'-','linewidth',lw,'color','blue');
hh2=plot(ta2,hosa2,'-','linewidth',lw,'color','red');
hh3=plot(ta3,hosa3,'-','linewidth',lw,'color',[0.9290, 0.6940, 0.1250]);
hh5=plot(tfo(1:250),hosfo(1:250),'-','linewidth',lw,'color','black');

% points=tvec+10;
% pointsy=.93*maxY;
% txt={'1','2','3','4','5','6'};
% text(10,pointsy,'PRE','fontsize',15);
% text(tvec(2)+5,pointsy,'LD','fontsize',15);
% for i=3:lt-1
%     text(points(i),pointsy,txt{i-2},'fontsize',15)
% end
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
set(gca,'ytick',[0:12000:96000]);
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
legend([hh0,hh4,hh1,hh2,hh3],'LDB','FO','B (12,000)','B (18,000)','B (24,000)','Position',[-0.245 0.325 1 1]);

pointsx=385;
pointsy=10^3*[4.5,13.5,19.5,25.5,82];
txt=[hlda,ha1,ha2,ha3,hfo];
for i=1:5
    text(pointsx,pointsy(i),['£' num2str(txt(i)) 'bn'],'fontsize',12);
end

hold off;