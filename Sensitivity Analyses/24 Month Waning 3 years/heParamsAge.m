function [phgs,pdgh,Threc,Thd]=heParamsAge(datax,ps)
%%

%nn=[3463,3726,3538,3260,3693,4022,4011,3926,3586,3919,4129,3890,3308,2982,2960,2069,1531+933+414+117+12];%England and Wales, mid-2019 estimate (ONS)
nn=datax.Npop';
nn=[nn(1:16),sum(nn(17:end))];%last age range in Knock et al. (2020) is 80+

ranges=[1,3,9,4];
nntot=[nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
nntot=repelem(nntot,ranges);
nnprop=nn./nntot;

subs=1:4;
subs=repelem(subs,ranges);

%%

% ps=     0.6;
ihr=    [0.030000 0.002600 0.000840	0.000420 0.000800 0.002600 0.004000	0.006300 0.012000 0.019000 0.023000	0.040000 0.096000 0.100000 0.240000	0.500000 0.500000];
phgs=   ihr./ps;

picu=   [0.2430 0.2890 0.3380 0.3890 0.4430 0.5030 0.5700 0.6530 0.7560 0.8660 0.9540 1.0000 0.9720 0.8540 0.6450 0.4020 0.1070];
pg=     1-picu;
pdicu=  [0.2820 0.2860 0.2910 0.2990 0.3100 0.3280 0.3530 0.3900 0.4460 0.5200 0.6040 0.7050 0.8060 0.8990 0.9690 1.0000 0.9180];
psd=    1-pdicu;

ifr=    [0.000310 0.000030 0.000010	0.000000 0.000010 0.000040 0.000060	0.000130 0.000310 0.000700 0.001160	0.002760 0.008670 0.012150 0.035120	0.084300 0.096960];
pdgh=   ifr./ihr;

phgs=   accumarray(subs',phgs.*nnprop);
pdgh=   accumarray(subs',pdgh.*nnprop);

%%

Tgr=    10.7;

Ttr=    2.5; 
Ticur=  15.6; 
Tsdr=   12.2;

Tgd=    10.3;

%Ttr=    2.5; 
Ticud=  11.8;

%Ttr=    2.5; 
Ticusd= 7; 
Tsdd=   8.1;

Threc=  Tgr*(pg*nn'/sum(nn))+(Ttr+Ticur+Tsdr)*(picu*nn'/sum(nn));
Thd=    Tgd*(pg*nn'/sum(nn))+(Ttr+Ticud)*((picu.*pdicu)*nn'/sum(nn))+(Ttr+Ticusd+Tsdd)*((picu.*psd)*nn'/sum(nn));

end