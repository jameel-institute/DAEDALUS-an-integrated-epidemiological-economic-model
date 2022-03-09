function f=heMakeDs(NN,x,data,int)

%NN - vector of populations including non-working.
%x - proportion of each sector open - not including non-working.
%The contact matrices are entirely determined by 'economics', both sector closures and working-from-home
%However, the succeptability of children may be included

%% COMMUNITY-COMMUNITY MATRIX:

childrenLessSus=0;

% if int==1
%     data.schoolA2=data.schoolA2*1;
% end
% if int==0
%     background=data.comm(1);
% elseif int==1
%     background=data.comm(1);
% else
%     background=data.comm(1);
% end

C=data.CM;%source: polymod-uk-all
Npop=data.Npop;
Npop(16)=sum(Npop(16:end));
Npop=Npop(1:16);

C=[C(:,1),sum(C(:,2:4),2),sum(C(:,5:13),2),sum(C(:,14:16),2)];%sum of the columns
C=[C(1,:);
   Npop(2:4)'*C(2:4,:)/sum(Npop(2:4));
   Npop(5:13)'*C(5:13,:)/sum(Npop(5:13));
   Npop(14:16)'*C(14:16,:)/sum(Npop(14:16))];%weighted average of the rows

Cav=([Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:16))]/sum(Npop))*sum(C,2);
C=data.comm*(C/Cav);

%%

lc=length(C);
lx=length(x);%Number of sectors%make this standard
ln=length(NN);

adInd=3;%Adult index
CworkRow=C(adInd,:);

%if length(NN)==lc+1%*1Dage
    %workContacts=27.4273;
    %f=[Cwork+x*workContacts,0;0,Cwork];
%else

NNrep=repmat(NN'/sum(NN),ln,1);%population proportion matrix
NNrel=NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));%proportion of adult population in each sector as vector

%Make A:
matA=zeros(ln,ln);
matA(lx+1:end,lx+1:end)=C;
matA(1:lx,lx+1:end)=repmat(CworkRow,lx,1);
matA(:,[1:lx,lx+adInd])=repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);

%% Modify depending on x:

if lx==1
    %Education:
    matA(lx+1,lx+1)=matA(lx+1,lx+1)+data.propschools*data.schoolA1*x(1);
    matA(lx+2,lx+2)=matA(lx+2,lx+2)+data.propschools*data.schoolA2*x(1);
    
    %Hospitality:
    %matA([1:lx,lx+3:ln],:)=matA([1:lx,lx+3:ln],:)+NNrep([1:lx,lx+3:ln],:)*datax.prophosp*datax.hospA34;
    matA([1:lx,lx+3],:)=matA([1:lx,lx+3],:)+NNrep([1:lx,lx+3],:)*data.prophosp*data.hospA3;
    matA(ln,:)=matA(ln,:)+NNrep(ln,:)*data.prophosp*data.hospA4;
    matA(lx+2,:)=matA(lx+2,:)+NNrep(lx+2,:)*data.prophosp*data.hospA2;%(1)
    
elseif lx==10
    %Education:
    matA(lx+1,lx+1)=matA(lx+1,lx+1)+data.propschools*data.schoolA1*x(9);
    matA(lx+2,lx+2)=matA(lx+2,lx+2)+data.propschools*data.schoolA2*x(9);
    
    %Hospitality:
    %matA([1:lx,lx+3:ln],:)=matA([1:lx,lx+3:ln],:)+NNrep([1:lx,lx+3:ln],:)*datax.prophosp*x(10)*datax.hospA34(1);
    matA([1:lx,lx+3],:)=matA([1:lx,lx+3],:)+NNrep([1:lx,lx+3],:)*data.prophosp*x(10)*data.hospA3(1);
    matA(ln,:)=matA(ln,:)+NNrep(ln,:)*data.prophosp*x(10)*data.hospA4(1);
    matA(lx+2,:)=matA(lx+2,:)+NNrep(lx+2,:)*data.prophosp*x(10)*data.hospA2;%(1)
      
elseif lx==36
    %Education:
    matA(lx+1,lx+1)=    matA(lx+1,lx+1)+    x(33)*  data.schoolA1;%mixing within age groups only
    matA(lx+2,lx+2)=    matA(lx+2,lx+2)+    x(33)*  data.schoolA2;
    
    %Hospitality:
    sects=[25,35];
    psub=data.NNs(sects);
    psub=sum(psub.*x(sects))/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
    matA([1:lx,lx+adInd],:)=    matA([1:lx,lx+adInd],:)+    psub*   data.hospA3*    NNrep([1:lx,lx+adInd],:);%mixing between all age groups, including pre-school
    matA(lx+2,:)=               matA(lx+2,:)+               psub*   data.hospA2*    NNrep(lx+2,:);
    matA(ln,:)=                 matA(ln,:)+                 psub*   data.hospA4*    NNrep(ln,:);

elseif lx==63
    %Education:
    matA(lx+1,lx+1)=    matA(lx+1,lx+1)+    x(55)*  data.schoolA1;%mixing within age groups only
    matA(lx+2,lx+2)=    matA(lx+2,lx+2)+    x(55)*  data.schoolA2;
    
    %Hospitality:
    sects=[36,58,59,60,62];
    psub=data.NNs(sects);
    psub=sum(psub.*x(sects))/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
    matA([1:lx,lx+adInd],:)=    matA([1:lx,lx+adInd],:)+    psub*   data.hospA3*    NNrep([1:lx,lx+adInd],:);%mixing between all age groups, including pre-school
    matA(lx+2,:)=               matA(lx+2,:)+               psub*   data.hospA2*    NNrep(lx+2,:);
    matA(ln,:)=                 matA(ln,:)+                 psub*   data.hospA4*    NNrep(ln,:);

else
    error('Unknown economic configuration!');
    
end

%Transport:
if int==0
    matA(1:lx,1:lx)=    matA(1:lx,1:lx)+    repmat(x',lx,1).*   data.travelA3(1).*  NNrep(1:lx,1:lx);%mixing between employed groups only
else
    matA(1:lx,1:lx)=    matA(1:lx,1:lx)+    repmat(x',lx,1).*   data.travelA3(1).*  NNrep(1:lx,1:lx).*  repmat(1-data.wfhAv,lx,1).*repmat(1-data.wfhAv',1,lx);%home-working has a compound effect
end

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES:

%Make B and C:
valB=data.B;
valC=data.C;
if int>0
    valB=valB.*(1-data.wfhAv).*(1-data.wfhAv);%home-working has a compound effect
    valC=valC.*(1-data.wfhAv);
end
valB(lx+1:ln)=0;
valC(lx:1:ln)=0;
x(lx+1:ln)=0;
matB=diag(x.*valB');
matC=repmat(x.*valC',1,ln).*NNrep;

%%

D=matA+matB+matC;

%%

if childrenLessSus==1
    %D(end-3,:)=.5*D(end-3,:);
    %D(end-2,:)=9/15*D(end-2,:);
    D(end-3,:)=(1-0.5)*D(end-3,:);
    D(end-2,:)=(1-(0.5*11/15))*D(end-2,:);
end

f=D;

end