function [g,h]=heSingleSim(xoptim)%,inp)

    load('UK63.mat','data');

    numInt=length(xoptim)/length(data.G);

    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt);%,inp);
    pr.sw=0;%switching off

    if numInt==3
        tvec=[-61.3373  87.6062  245  306  367  426];
        %tvec=[-60.7399  86.8162  245  306  367  426];%ChLessSus fit
    elseif numInt==6
        %tvec=[-61.3373  87.6062  245  275  306  336  367  398  426];%6x1
        tvec=[-61.3373  87.6062  245  306  367  426  487  548  610];%6x2
    elseif numInt==12
        tvec=[-61.3373  87.6062  245  306  367  426  487  548  610  671  732  791  852  913  975];%12x2
    else
        error('Input not consistent!');
    end
    % tvec=[1  linspace(pr.Tm,  365*2+1,  numInt+1)];
    
    woptim=xoptim.^(1/pr.a);
    
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,woptim,tvec,0,data);

    plotMultiOutd_sch(f,xoptim,tvec,data);

    format bank;
    h(1)=f(end,5);%deaths
    % fullx=[ones(length(data.G),1),reshape(xoptim,length(data.G),numInt)];
    fullx=[ones(length(data.G),1),data.xmin',reshape(xoptim,length(data.G),numInt)];
    dgva=(12/365)*data.obj;
    h(2)=((tvec(end)-tvec(1))*sum(dgva)-(tvec(2:end)-tvec(1:end-1))*sum(fullx.*dgva,1)')/1000;%GDP loss ($, billion)
    
    h(2)=round(h(2),2);
    % figure(1);
    % title(['GDP loss \$' num2str(h(2)) ' (billion)']);
    
    %%% Schools Paper:
    h=sum(xoptim.*repmat((24/numInt)*data.obj,numInt,1))/1000;
    %%%

end