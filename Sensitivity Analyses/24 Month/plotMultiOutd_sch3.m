load('xoptim.mat')
[g,h]=heSingleSim(xoptim)

close figure 2;
figure(1);
hold on;
legend('AutoUpdate','off')

load('../../R1/12 Month/xoptim.mat');
    load('UK63.mat','data');
    numInt=length(xoptim)/length(data.G);
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt);%,inp);
    pr.sw=0;%switching off
    tvec=[-61.3373   87.6062  245 306  367  426  487  548  610];%6x2
    woptim=xoptim.^(1/pr.a);
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,woptim,tvec,0,data);
    numPeriods=length(tvec)-1;
    numSectors=length(data.G);
    lt=length(tvec);
    dodiff=1;
    t1=f(:,1);  
    I1=f(:,3);
    h1=f(:,4);
    scal=sum(data.Npop)/(10*10^6);
    lw=2.5;
    hh2=plot(t1(t1>tvec(3)),I1(t1>tvec(3))/scal,'-.','linewidth',lw,'color',[1.00,0.50,0]);
    hh3=plot(t1(t1>tvec(3)),h1(t1>tvec(3)),'-.','linewidth',lw,'color',[0.38,0,0.50]);
    
load('../../Original/B2.mat');
    numInt=length(xoptim)/length(data.G);
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt);%,inp);
    pr.sw=0;%switching off
    tvec=[-61.3373   87.6062  245 306  367  426];    woptim=xoptim.^(1/pr.a);
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,woptim,tvec,0,data);
    numPeriods=length(tvec)-1;
    numSectors=length(data.G);
    lt=length(tvec);
    dodiff=1;
    t1=f(:,1);  
    I1=f(:,3);
    h1=f(:,4);
    scal=sum(data.Npop)/(10*10^6);
    lw=2.5;
    hh2=plot(t1(t1>tvec(3)),I1(t1>tvec(3))/scal,'--','linewidth',lw,'color',[1.00,0.50,0]);
    hh3=plot(t1(t1>tvec(3)),h1(t1>tvec(3)),'--','linewidth',lw,'color',[0.38,0,0.50]);