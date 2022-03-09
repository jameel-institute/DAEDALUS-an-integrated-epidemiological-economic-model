close all;
clear all;

addpath('../');
addpath('../Output Data/GDP-63/');

load('UK63.mat','data');
load('A2.mat');
numInt=length(xoptim)/length(data.G);
tvec=[-61.3373   87.6062  245 306  367  426];

wfh=data.wfhAv;

for i=1:11
    
    data.wfhAv=wfh*(1-(i-1)*0.02);
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt);%,inp);
    pr.sw=0;%switching off

    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
    
    h(:,2*i-1)= f(:,1);
    h(:,2*i)=   f(:,4);

end

%%

numPeriods=length(tvec)-1;
numSectors=length(data.G);
lt=length(tvec);
dodiff=1;

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
lw=2.0;
cmap=lines(2);
thresh=[12000,18000,24000];
numThresh=length(thresh);

maxY=48000;
%
hold on;
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=2
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end

cmap=jet(11);

for i=11:-1:1
    
    plot(h(:,2*i-1),h(:,2*i),'-','linewidth',lw,'color',cmap(i,:));

end

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

% pointsx=385;
% pointsy=10^3*[2.5,13.5,19.5,25.5,82];
% txt=[hlda,ha1,ha2,ha3,hfo];
% for i=1:5
%     text(pointsx,pointsy(i),['£' num2str(txt(i)) 'bn'],'fontsize',12);
% end

hold off;

%%

load('B2.mat');

for i=1:11
    
    data.wfhAv=wfh*(1-(i-1)*0.02);
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt);%,inp);
    pr.sw=0;%switching off

    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);
    
    if i~=11
        hh(:,2*i-1)= f(:,1);
        hh(:,2*i)=   f(:,4);
    else
        hhh(:,1)=    f(:,1);   
        hhh(:,2)=    f(:,4);
    end
    
end

%%

numPeriods=length(tvec)-1;
numSectors=length(data.G);
lt=length(tvec);
dodiff=1;

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
lw=2.0;
cmap=lines(2);
thresh=[12000,18000,24000];
numThresh=length(thresh);

maxY=48000;
%
hold on;
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=2
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end

cmap=jet(11);

plot(hhh(:,1),hhh(:,2),'-','linewidth',lw,'color',cmap(11,:));
for i=10:-1:1
    
    plot(hh(:,2*i-1),hh(:,2*i),'-','linewidth',lw,'color',cmap(i,:));

end

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

% pointsx=385;
% pointsy=10^3*[2.5,13.5,19.5,25.5,82];
% txt=[hlda,ha1,ha2,ha3,hfo];
% for i=1:5
%     text(pointsx,pointsy(i),['£' num2str(txt(i)) 'bn'],'fontsize',12);
% end

hold off;