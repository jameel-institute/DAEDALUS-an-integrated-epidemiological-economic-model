function f=plotMultiOutd_sch2(x1,x2,tvec,data)
numPeriods=length(tvec)-1;
numSectors=length(data.G);
lt=length(tvec);
dodiff=1;

if numPeriods==5
	xlabels2=({'PRE','LD','1','2','3'});
elseif numPeriods==8
    %xlabels2=({'Jan','Mar 26th','Sep','Nov','Jan'});
    xlabels2=({'PRE','LD','1','2','3','4','5','6'});
elseif numPeriods==14
    %xlabels2=({'Jan','Mar 26th','Sep','Nov','Jan'});
    xlabels2=({'Jan','LD','Sep','Nov','Jan','Mar','May','Jul','Sep','Nov','Jan','Mar','May','Jul'});
else
    error('Data missing for nunmPeriods')
end

f=figure('Units','centimeters','Position',[0 0 20 18]);
fs=15;set(f,'DefaultAxesFontSize',15);
x1=[ones(numSectors,1);data.xmin';x1];
x2=[ones(numSectors,1);data.xmin';x2];
x=reshape(x1-x2,numSectors,numPeriods);
xvec=(1:numPeriods)';
xvec=xvec-.5;
yvec=(1:10:numSectors)';
ylab=num2str(yvec);
yvec=yvec-.5;
%colormap gray
b=imagesc(x);
set(b,'AlphaData',~isnan(x))
set(gca,'YDir','normal')
xlabel('Time','FontSize',fs)
ylabel('Economic Sector','FontSize',fs);
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.2 0 0]);
set(gca,'xtick',xvec,'xticklabels',xlabels2,'ytick',yvec,'yticklabels',ylab);%{'PRE','LD',xlab(3:end,:)}
axis square;%([0,numPeriods,.5,63.5])%yvec(1),yvec(end)+1])
xtickangle(45);
colormap(flipud(customcolormap_preset('pink-white-green')));
caxis([-1,1])
colorbar
grid on
grid minor
ax = gca;
ax.XAxis.MinorTickValues = 1.5:1:5.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 1.5:1:63.5;
box on

end