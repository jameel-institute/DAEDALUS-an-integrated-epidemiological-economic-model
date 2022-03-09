function f=plotMultiOutd_sch(f1,x1,tvec,data)
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
t1=f1(:,1); 
s1=f1(:,2); 
I1=f1(:,3);
h1=f1(:,4);
d1=f1(:,5);
v1=f1(:,6);%heSimCovid19 output
if dodiff==1
    inc1=-diff(s1,1);%f1=Sout
    tdiff=diff(t1,1);
    inc1=inc1./tdiff;%repmat - older version?
end

maxY=max([inc1;I1]);%inc2;h2
maxY=48000;
%
hold on
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end
for j=1:numThresh
    plot([0,tvec(end)],[thresh(j),thresh(j)],'-','linewidth',lw,'color',.5*[1,1,1])
end

scal=sum(data.Npop)/(10*10^6);
% hh1=plot(t1(1:end-1)+0.5*tdiff,inc1,'-','linewidth',lw,'color','yellow');
hh2=plot(t1,I1/scal,'-','linewidth',lw,'color',[1.00,0.50,0]);
hh3=plot(t1,h1,'-','linewidth',lw,'color',[0.38,0,0.50]);
% hh4=plot(t1,d1,'-','linewidth',lw,'color','black');
% hh5=plot(t1,v1,'-','linewidth',lw,'color','green');

%}
points=tvec+10;
pointsy=.93*maxY;
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
ylabel('Number','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-20 0 0]);
%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[1,32,61,92,122,153,183,214,245,275,306,336,367,398]);
set(gca,'ytick',[0:6000:48000]);
set(gca,'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'});
if numPeriods==5
	xlabels2=({'Jan','LD','Sep','Nov','Jan'});
elseif numPeriods==8
    %xlabels2=({'Jan','Mar 26th','Sep','Nov','Jan'});
    xlabels2=({'Jan','LD','Sep','Oct','Nov','Dec','Jan','Feb'});
else
    error('Data missing for nunmPeriods')
end
xtickangle(45);
ax = gca;
ax.YAxis.Exponent = 3;
box on;
grid on;
grid minor;
%legend([hh2,hh3],'Prevalence (per 10m)','Hospital Occupancy','Position',[0.2275 0.273 1 1]);
legend([hh2,hh3],'Prevalence (per 10m)','Hospital Occupancy','Position',[0.2275 0.380 1 1]);
%legend([hh1,hh2,hh3,hh4,hh5],'Incidence','I','H','D','V','location','northwest');
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
hold off;

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',15);
x1=[ones(numSectors,1);data.xmin';x1];
x=reshape(x1,numSectors,numPeriods);
xvec=(1:numPeriods)';
xvec=xvec-.5;
yvec=(1:10:numSectors)';
ylab=num2str(yvec);
yvec=yvec-.5;
%colormap gray
imagesc(x)
set(gca,'YDir','normal')
xlabel('Time','FontSize',fs)
ylabel('Economic Sector','FontSize',fs);
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.2 0 0]);
set(gca,'xtick',xvec,'xticklabels',xlabels2,'ytick',yvec,'yticklabels',ylab);%{'PRE','LD',xlab(3:end,:)}
axis square;%([0,numPeriods,.5,63.5])%yvec(1),yvec(end)+1])
xtickangle(45);
%colormap(autumn);
caxis([0,1])
colorbar
grid on
grid minor
ax = gca;
ax.XAxis.MinorTickValues = 1.5:1:5.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 1.5:1:63.5;
box on

%%

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
lw=2.5;
cmap=lines(67);
t1=f1(:,1); 
IS=f1(:,6+1:6+67);
ISt=f1(:,end-2);

maxY=1;
%
hold on
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end

for i=1:63
    plot(t1,100*IS(:,i),'-','linewidth',0.5,'color',cmap(i,:));
end
hh2=plot(t1,100*ISt,'-','linewidth',lw,'color','black');
hh3=plot(t1,100*IS(:,67),'-.','linewidth',lw,'color','black');

%}
% points=tvec+10;
% pointsy=.85*maxY;
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
ylabel('Percentage of Sector (%)','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-20 0 0]);
%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[1,32,61,92,122,153,183,214,245,275,306,336,367,398]);
set(gca,'ytick',[0:0.1:100]);
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
%ax.YAxis.Exponent = 1;
box on;
grid on;
grid minor;
legend([hh2,hh3],'Working Age','Retired Age','location','northeast');
%legend([hh1,hh2,hh3,hh4,hh5],'Incidence','I','H','D','V','location','northwest');
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
hold off;
title('Symptomatic Infections');

%%

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
lw=2.5;
cmap=lines(67);
t1=f1(:,1); 
HH=f1(:,6+67+1:6+67+67);
HHt=f1(:,end-1);

maxY=1;
%
hold on
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end

for i=1:63
    plot(t1,100*HH(:,i),'-','linewidth',0.5,'color',cmap(i,:));
end
hh2=plot(t1,100*HHt,'-','linewidth',lw,'color','black');
hh3=plot(t1,100*HH(:,67),'-.','linewidth',lw,'color','black');

%}
% points=tvec+10;
% pointsy=.85*maxY;
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
ylabel('Percentage of Sector (%)','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-20 0 0]);
%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[1,32,61,92,122,153,183,214,245,275,306,336,367,398]);
set(gca,'ytick',[0:0.1:100]);
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
%ax.YAxis.Exponent = 1;
box on;
grid on;
grid minor;
legend([hh2,hh3],'Working Age','Retired Age','location','northeast');
%legend([hh1,hh2,hh3,hh4,hh5],'Incidence','I','H','D','V','location','northwest');
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
hold off;
title('Hospitalised Infections');

%%

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',fs);
lw=2.5;
cmap=lines(67);
t1=f1(:,1); 
DD=f1(:,6+2*67+1:6+2*67+67);
DDt=f1(:,end);

maxY=1;
%
hold on
for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
end

for i=1:63
    plot(t1,100*DD(:,i),'-','linewidth',0.5,'color',cmap(i,:));
end
hh2=plot(t1,100*DDt,'-','linewidth',lw,'color','black');
hh3=plot(t1,100*DD(:,67),'-.','linewidth',lw,'color','black');

%}
% points=tvec+10;
% pointsy=.85*maxY;
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
ylabel('Percentage of Sector (%)','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-20 0 0]);
%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[1,32,61,92,122,153,183,214,245,275,306,336,367,398]);
set(gca,'ytick',[0:0.1:100]);
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
%ax.YAxis.Exponent = 1;
box on;
grid on;
grid minor;
legend([hh2,hh3],'Working Age','Retired Age','location','northeast');
%legend([hh1,hh2,hh3,hh4,hh5],'Incidence','I','H','D','V','location','northwest');
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
hold off;
title('Dead');

end