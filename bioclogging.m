function biologicalclogging

load Sn;
tt=logspace(-2,6,100)';     % Temporal step,[s]
Q=50/86400;                % Injection rate,[m3/s]
m=6;                       % Aquifer thickness,[m]
re=1000;                    % Distance that is far enough,[m]
rw=0.001;                  % Well radius, [m]
rc=0.001;                  % Well casing radius,[m]

N=100;                    % Spatial node numbers
poroini=0.33;                  % Porosity,[-]
kini=1e-5;                    % Initial hydraulic conductivity, [m/s]
alpha=0.1;
Ss=0.001;                    % Storage coefficient,[-]
CO_inj=8;                   % Injection DO concentration, [g/m3]
Csub_inj=0;
Cbio_inj=0;
CeO=0.01;                   % Oxygen concentration far away from well,[g/m3]
Cesub=0;
Cebio=0;

Csub=100;                  % Initial substrate concentration,[g/m3]
Csub_B=Csub*0.05/100;
Csub_T=Csub*6.6/100;
Csub_X=Csub*6.9/100;
rho_bulk=1.77e6;            % Bulk density, [g/m3]
rhob=1000;                 % Density of biomass [g/m3 dry mass/wet volume] Newcomer et al. (2016)
katt=3.2778e-4;               % Attatchment rate [1/s] Farrokhian Firouzi et al., 2015
kdet=8.9444e-6;               % Detachment rate [1/s] Farrokhian Firouzi et al., 2015
C_bio_m=0.2;                   % Initial mobile biomass concentration, [g/m3 dry mass/wet volume]
C_bio_im=C_bio_m*katt/kdet/poroini/rho_bulk;  % Initial immobile biomass concentration (dry mass), [g/m3]
L_cont=10;
b=0;                 % Decay rate, [s-1] Kim et al. [23]
c=0.8;                     % Thullner et al. (2002b)
kmin=0.008;                % Newcomer et al. (2016)
kNAPL=5e-3/86400;          %Mass-transfer coefficient of NAPL to aqueous phase [Brauner and Widdowson, 2001]

pars.kNAPL=kNAPL;
pars.poroini=poroini;
pars.Ss=Ss;
pars.alpha=alpha;
pars.Q=Q;
pars.m=m;
pars.re=re;
pars.rw=rw;
pars.rc=rc;
pars.CO_inj=CO_inj;
pars.Csub_inj=Csub_inj;
pars.Cbio_inj=Cbio_inj;
pars.rho_bulk=rho_bulk;
pars.rhob=rhob;
pars.katt=katt;
pars.kdet=kdet;
pars.b=b;
pars.CeO=CeO;
pars.Cesub=Cesub;
pars.Cebio=Cebio;
%***************** bioreaction parameters ************************************
KO=0.2;                   % Half-saturation constant of oxygen, [g/m3]
b=0.23/86400;               % Biomass specific loss rate, [1/s]
De=1e-9;                   % Molecular diffusion rate of substrate through water,[m2/s]
pars.KO=KO;
pars.b=b;
pars.De=De;

%**************************************************************************

%Develop r-axis
i=[0:N]';
rBND=10.^(log10(rw)+i*(log10(re)-log10(rw))/N);   % Spatial step, [m]
r=(rBND(1:N)+rBND(2:N+1))/2; 
rIn=rBND(1:N);
rOut=rBND(2:N+1);
pars.rIn=rIn;
pars.rOut=rOut;
pars.r=r;
pars.rBND=rBND;
loc_cont=round((log10(L_cont)-log10(rw))/((log10(re)-log10(rw))/N));
%Define initial condition
CO=CeO*ones(N,1);         % Fixed concentration in well, so "N-1" is defined
Csub_B=[Csub_B*ones(loc_cont,1);zeros(N-loc_cont,1)];
Csub_T=[Csub_T*ones(loc_cont,1);zeros(N-loc_cont,1)];
Csub_X=[Csub_X*ones(loc_cont,1);zeros(N-loc_cont,1)];
Cbio_m=[C_bio_m*ones(loc_cont,1);zeros(N-loc_cont,1)];
Cbio_im=[C_bio_im*ones(loc_cont,1);zeros(N-loc_cont,1)];
N=N+1;                     % The "1" is defined for well
pars.N=N;
s0=zeros(N,1);             % Initial drawdown, [m]

h0_iTime_old=[s0;CO;Csub_B;Csub_T;Csub_X;Cbio_m;Cbio_im];

t1=tt;
%********************** Using ODE under FOR LOOP   ************************
%********** Size setup for each parameter****************
Time_num=numel(t1);

poro=zeros(Time_num,N-1);
k=zeros(Time_num,N-1);

nbio=Cbio_im'*rho_bulk/rhob;
poroini=poroini*ones(1,N-1);
poro(1,:)=poroini-nbio;

poro_rel=(poroini-nbio)./poroini;
kn=(((poro_rel-poroini)./(1-poroini)).^c+kmin)*(1/(1+kmin));
kini=kini*ones(1,N-1);
k(1,:)=kn.*kini;

% Define empty matrix to store each parameter
ss=zeros(Time_num,N);
CO=zeros(Time_num,N-1);
CsubB=zeros(Time_num,N-1);
CsubT=zeros(Time_num,N-1);
CsubX=zeros(Time_num,N-1);
Cbio_m=zeros(Time_num,N-1);
Cbio_im=zeros(Time_num,N-1);
CsubTOT=zeros(Time_num,N-1);
%*********************************************************
for iTime=1:Time_num
    if iTime==1
        t=[0,t1(iTime)]';
        h0=h0_iTime_old;
        
    pars.poro=poro(iTime,:)';
    pars.k=k(iTime,:)';
    else
        t=[0,t1(iTime)-t1(iTime-1)]';
        h0=h0_iTime_old;
    
    pars.poro=poro(iTime,:)';
    pars.k=k(iTime,:)';
    end
    
options=[];
[tt,h]=ode15s(@odefun,t,h0,options,pars);
ss(iTime,:)=h(end,1:N);
CO(iTime,:)=h(end,N+1:2*N-1);
CsubB(iTime,:)=h(end,2*N:3*N-2);
CsubT(iTime,:)=h(end,3*N-1:4*N-3);
CsubX(iTime,:)=h(end,4*N-2:5*N-4);
Cbio_m(iTime,:)=h(end,5*N-3:6*N-5);
Cbio_im(iTime,:)=h(end,6*N-4:end);
CsubTOT(iTime,:)=CsubB(iTime,:)+CsubT(iTime,:)+CsubX(iTime,:);

% The drawdown and oxygen concentration results calculated from the last 
% step of ODE function is used as the next initial value of next ODE
h0_iTime_old=[ss(iTime,:),CO(iTime,:),CsubB(iTime,:),CsubT(iTime,:),CsubX(iTime,:),Cbio_m(iTime,:),Cbio_im(iTime,:)]';  

%********************** Calculating bioclogging   *************************
nbio=Cbio_im(iTime,:)*rho_bulk/rhob;
poro(iTime+1,:)=poroini-nbio;

poro_rel=poro(iTime+1,:)./poroini;
kn=(((poro_rel-poroini)./(1-poroini)).^c+kmin)*(1/(1+kmin));
k(iTime+1,:)=kn.*kini;
end
%*****************************************************************
hw=ss(:,1);
poro(1,:)=[];
k(1,:)=[];
% save hwbio hw -ascii;
save hwCsub2 hw -ascii;
dporo=1-poro/poroini(1);
dk=1-k/kini(1);

figure(4)
loglog(t1,hw,'r','Linewidth',2);
xlabel('\it t \rm[s]','FontName','Times New Roman','Fontweight','bold','FontSize',23);
ylabel('\it s \rm[m]','FontName','Times New Roman','Fontweight','bold','FontSize',23);
xlim([t1(1) t1(end) ]);
% ylim([1e-3 20]);
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',15);
title('Drawdown','Fontsize',12);
hold on;

figure(2)
[x,y]=meshgrid(r,t1);
surf(x,y,dporo);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(get(gca,'XLabel'),'Rotation',12.5);
set(get(gca,'YLabel'),'Rotation',-18);
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',16);
xlabel('Distance [\itm\rm]','FontName','Times New Roman','Fontweight','bold','FontSize',18,'position',[0.2 0.0002]);
ylabel('Time [\its\rm]','FontName','Times New Roman','Fontweight','bold','FontSize',18,'position',[1e-4 5]);
zlabel('Reduction of porosity [-]','FontName','Times New Roman','Fontweight','bold','FontSize',18);
set(gca,'xtick',[0.001 0.01 0.1 1 10 100 1000]);
set(gca,'ytick',[1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5 1e6]);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','points');
set(gcf,'PaperPosition',[0 0 1040 480]);
colorbar('position',[0.94 0.1 0.03 0.8]);
% set(gca,'XMinorGrid','on');
% set(gca,'YMinorGrid','on');
xlim([rc re]);
ylim([t1(1) t1(end)]);
zlim([0 1]);
shading flat;
print('-djpeg100','-r700','poro');

figure(3)
[x,y]=meshgrid(r,t1);
surf(x,y,dk);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(get(gca,'XLabel'),'Rotation',12.5);
set(get(gca,'YLabel'),'Rotation',-18);
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',16);
xlabel('Distance [\itm\rm]','FontName','Times New Roman','Fontweight','bold','FontSize',18,'position',[0.2 0.0002]);
ylabel('Time [\its\rm]','FontName','Times New Roman','Fontweight','bold','FontSize',18,'position',[1e-4 5]);
zlabel('Reduction of hydraulic conductivity [-]','FontName','Times New Roman','Fontweight','bold','FontSize',18);
set(gca,'xtick',[0.001 0.01 0.1 1 10 100 1000]);
set(gca,'ytick',[1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5 1e6]);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','points');
set(gcf,'PaperPosition',[0 0 1040 480]);
colorbar('position',[0.94 0.1 0.03 0.8]);
% set(gca,'XMinorGrid','on');
% set(gca,'YMinorGrid','on');
xlim([rc re]);
ylim([t1(1) t1(end)]);
zlim([0 1]);
shading flat;
% print('-djpeg100','-r700','k');



%**************************************************************************
%**************************************************************************
% function dhdt=odefun(t,h,pars,wait,tMax)
function dhdt=odefun(t,h,pars)
N=pars.N;
s(1:N,:)=h(1:N,:);
CO(1:N-1,:)=h(N+1:2*N-1,:);
CsubB(1:N-1,:)=h(2*N:3*N-2,:);
CsubT(1:N-1,:)=h(3*N-1:4*N-3,:);
CsubX(1:N-1,:)=h(4*N-2:5*N-4,:);
Cbio_m(1:N-1,:)=h(5*N-3:6*N-5,:);
Cbio_im(1:N-1,:)=h(6*N-4:end,:);

sw=s(1,:);
s(1,:)=[];

Q=pars.Q;
m=pars.m;
r=pars.r;
rIn=pars.rIn;
rOut=pars.rOut;
re=pars.re;
rw=pars.rw;
rc=pars.rc;
k=pars.k;
Ss=pars.Ss;
alpha=pars.alpha;
poro=pars.poro;
CeO=pars.CeO;
Cesub=pars.Cesub;
Cebio=pars.Cebio;
CO_inj=pars.CO_inj;
Csub_inj=pars.Csub_inj;
Cbio_inj=pars.Cbio_inj;
rho_bulk=pars.rho_bulk;
kNAPL=pars.kNAPL;

N=N-1;
%=======================    flow field   =================================
%Hydraulic gradient in r-dir (general condition)
i=2:N;

dsdrIn(i,:)=(s(i-1,:)-s(i,:))./(r(i)-r(i-1));
dsdrOut(i-1,:)=dsdrIn(i,:);

%Apply well-head condition to near-side boundary
dsdrIn(1,:)=(sw-s(1,:))./(r(1)-rw);

%Apply fixed-head condition to far-side boundary
dsdrOut(N,:)=(s(N,:)-0)./(re-r(N));

%Caculate water velocity using Darcy's Law
qIn=k.*sign(dsdrIn).*abs(dsdrIn);
qOut=k.*sign(dsdrOut).*abs(dsdrOut);

%Apply continuity in aquifer
% dsdt=(rIn.*qIn-rOut.*qOut)./r./(rOut-rIn)*m/Ss;
dsdt=(rIn.*qIn-rOut.*qOut)./r./(rOut-rIn)*m./Ss;

%Apply well-bore equation
dswdt=(Q-2*pi*rw*m*qIn(1,:))/(pi*rc^2);

dsdt=[dswdt;dsdt];
%=======================  Solute transport  ==============================

%           ============  Advection  ===============
% OXYGEN
dCdtO_adv(i,:)=qIn(i,:)./poro(i,:).*(CO(i-1,:)-CO(i,:))./(r(i)-r(i-1));
dCdtO_adv(1,:)=qIn(1,:)./poro(1,:).*(CO_inj-CO(1,:))./(r(1)-rw);
% Benzene
dCdtB_adv(i,:)=qIn(i,:)./poro(i,:).*(CsubB(i-1,:)-CsubB(i,:))./(r(i)-r(i-1));
dCdtB_adv(1,:)=qIn(1,:)./poro(1,:).*(Csub_inj-CsubB(1,:))./(r(1)-rw);
% Toluene
dCdtT_adv(i,:)=qIn(i,:)./poro(i,:).*(CsubT(i-1,:)-CsubT(i,:))./(r(i)-r(i-1));
dCdtT_adv(1,:)=qIn(1,:)./poro(1,:).*(Csub_inj-CsubT(1,:))./(r(1)-rw);
% Xylene
dCdtX_adv(i,:)=qIn(i,:)./poro(i,:).*(CsubX(i-1,:)-CsubX(i,:))./(r(i)-r(i-1));
dCdtX_adv(1,:)=qIn(1,:)./poro(1,:).*(Csub_inj-CsubX(1,:))./(r(1)-rw);
% BIOMASS
dCdtBIO_adv(i,:)=qIn(i,:)./poro(i,:).*(Cbio_m(i-1,:)-Cbio_m(i,:))./(r(i)-r(i-1));
dCdtBIO_adv(1,:)=qIn(1,:)./poro(1,:).*(Cbio_inj-Cbio_m(1,:))./(r(1)-rw);

%           ============  Dispersion  ===============
%Calculate influx at each inter-face
% OXYGEN
dCdrOIn_dis(i,:)=alpha*abs(qIn(i,:))./poro(i,:).*(CO(i-1,:)-CO(i,:))./(r(i)-r(i-1));
dCdrOOut_dis(i-1,:)=dCdrOIn_dis(i,:);
% Benzene
dCdrBIn_dis(i,:)=alpha*abs(qIn(i,:))./poro(i,:).*(CsubB(i-1,:)-CsubB(i,:))./(r(i)-r(i-1));
dCdrBOut_dis(i-1,:)=dCdrBIn_dis(i,:);
% Toluene
dCdrTIn_dis(i,:)=alpha*abs(qIn(i,:))./poro(i,:).*(CsubT(i-1,:)-CsubT(i,:))./(r(i)-r(i-1));
dCdrTOut_dis(i-1,:)=dCdrTIn_dis(i,:);
% Xylene
dCdrXIn_dis(i,:)=alpha*abs(qIn(i,:))./poro(i,:).*(CsubX(i-1,:)-CsubX(i,:))./(r(i)-r(i-1));
dCdrXOut_dis(i-1,:)=dCdrXIn_dis(i,:);
% BIOMASS
dCdrBIOIn_dis(i,:)=alpha*abs(qIn(i,:))./poro(i,:).*(Cbio_m(i-1,:)-Cbio_m(i,:))./(r(i)-r(i-1));
dCdrBIOOut_dis(i-1,:)=dCdrBIOIn_dis(i,:);

%Calculate flux at each boundary
% OXYGEN
dCdrOIn_dis(1,:)=0;
dCdrOOut_dis(N,:)=alpha*abs(qOut(N,:))./poro(N,:).*(CO(N,:)-CeO)./(re-r(N));
dCdtO_dis=(rIn.*dCdrOIn_dis-rOut.*dCdrOOut_dis)./(rOut-rIn)./r;
% Benzene
dCdrBIn_dis(1,:)=0;
dCdrBOut_dis(N,:)=alpha*abs(qOut(N,:))./poro(N,:).*(CsubB(N,:)-Cesub)./(re-r(N));
dCdtB_dis=(rIn.*dCdrBIn_dis-rOut.*dCdrBOut_dis)./(rOut-rIn)./r;
% Toluene
dCdrTIn_dis(1,:)=0;
dCdrTOut_dis(N,:)=alpha*abs(qOut(N,:))./poro(N,:).*(CsubT(N,:)-Cesub)./(re-r(N));
dCdtT_dis=(rIn.*dCdrTIn_dis-rOut.*dCdrTOut_dis)./(rOut-rIn)./r;
% Xylene
dCdrXIn_dis(1,:)=0;
dCdrXOut_dis(N,:)=alpha*abs(qOut(N,:))./poro(N,:).*(CsubX(N,:)-Cesub)./(re-r(N));
dCdtX_dis=(rIn.*dCdrXIn_dis-rOut.*dCdrXOut_dis)./(rOut-rIn)./r;
% BIOMASS
dCdrBIOIn_dis(1,:)=0;
dCdrBIOOut_dis(N,:)=alpha*abs(qOut(N,:))./poro(N,:).*(Cbio_m(N,:)-Cebio)./(re-r(N));
dCdtBIO_dis=(rIn.*dCdrBIOIn_dis-rOut.*dCdrBIOOut_dis)./(rOut-rIn)./r;

%========================bioreaction and source sink=======================
rhob=pars.rhob;
b=pars.b;
katt=pars.katt;
kdet=pars.kdet;
poroini=pars.poroini;

% Maximum organic utilization rate
miu_B=8.4/86400;        % Maximum utilization rate of Benzene [1/s] Freitas et al., 2011; Molson et al., 2002
miu_T=10.68/86400;      % Maximum utilization rate of Toluene [1/s] Freitas et al., 2011; Bekins et al., 1998
miu_X=6/86400;          % Maximum utilization rate of Xylene [1/s] Freitas et al., 2011; Bekins et al., 1998
%miu_B=8.4/86400;
%miu_T=10.68/86400;
%miu_X=6/86400;   
% Microbial yield coefficient
Y_B=1.5;                % Microbial yield coefficient of Benzene [mg biomass/mg substrate] Bekins et al., 1998
Y_T=1.22;               % Microbial yield coefficient of Toluene [mg biomass/mg substrate] Bekins et al., 1998
Y_X=1.3;                % Microbial yield coefficient of Xylene [mg biomass/mg substrate]  Bekins et al., 1998
% Organic half rate utilization concentration
Ks_B=0.3;               % Organic half rate utilization concentration of Benzene [g/m3] Bekins et al., 1998
Ks_T=0.1;               % Organic half rate utilization concentration of Benzene [g/m3] Bekins et al., 1998 
Ks_X=0.0007;            % Organic half rate utilization concentration of Xylene [g/m3] Bekins et al., 1998 
% Oxygen half rate utilization concentration
KO=2;                   % Freitas et al., 2011; Molson et al., 2002         
% Ratio of oxygen to substrate consumed
F_B=3.08;               % Benzene [mg Oxygen/mg sub]
F_T=3.13;               % Toluene [mg Oxygen/mg sub]
F_X=3.17;               % Xylene [mg Oxygen/mg sub]
% Partition coefficient
Kd_B=1.3155e-5*1.01;            % Distribution coefficient of Benzene, [m3/g] "Ma and Li, 2014"Boggs et al. (1993); Brauner and Widdowson, 2001
Kd_T=1e-5*1.01;            % Distribution coefficient of Toluene, [m3/g] Boggs et al. (1993); Brauner and Widdowson, 2001
Kd_X=8.5e-6*1.01;            % Distribution coefficient of p-Xylene, [m3/g] Boggs et al. (1993); Brauner and Widdowson, 2001
%Kd_B=1.3155e-5;
%Kd_T=1e-5;
%Kd_X=8.5e-6;
% Equilibrium concentration of each PHC
Csol_B=1790*0.000834891;
Csol_T=500*0.093453407;
Csol_X=185*0.084332603;
% Retardation factor
Rd_B=1+rho_bulk./poro*Kd_B;             % Retardation factor of Benzene
Rd_T=1+rho_bulk./poro*Kd_T;             % Retardation factor of Toluene
Rd_X=1+rho_bulk./poro*Kd_X;             % Retardation factor of p-Xylene

% bioreaction rate
limi=1-Cbio_im*rho_bulk/rhob/poroini;  % Kildsgaard and Engesgaard, 2001

rBm=miu_B*Cbio_m.*(CsubB./(CsubB+Ks_B)).*(CO./(CO+KO)).*limi;
rTm=miu_T*Cbio_m.*(CsubT./(CsubT+Ks_T)).*(CO./(CO+KO)).*limi;
rXm=miu_X*Cbio_m.*(CsubX./(CsubX+Ks_X)).*(CO./(CO+KO)).*limi;

rBim=miu_B*Cbio_im.*(CsubB./(CsubB+Ks_B)).*(CO./(CO+KO)).*limi;
rTim=miu_T*Cbio_im.*(CsubT./(CsubT+Ks_T)).*(CO./(CO+KO)).*limi;
rXim=miu_X*Cbio_im.*(CsubX./(CsubX+Ks_X)).*(CO./(CO+KO)).*limi;

% Biomass, substrate and DO variation rates
rBIO_m=rBm+rTm+rXm;
rBIO_im=rBim+rTim+rXim;
r_B=-(rBm+rBim*rho_bulk./poro)/Y_B;
r_T=-(rTm+rTim*rho_bulk./poro)/Y_T;
r_X=-(rXm+rXim*rho_bulk./poro)/Y_X;
rO=r_B*F_B+r_T*F_T+r_X*F_X;

%Calculate concentration variation due to bioreaction-advection-dispersion
% OXYGEN
dCdtO=dCdtO_adv+dCdtO_dis+rO;
% Benzene
dCdt_B=(1./Rd_B).*(dCdtB_adv+dCdtB_dis+r_B+kNAPL*(Csol_B-CsubB));
% Toluene
dCdt_T=(1./Rd_T).*(dCdtT_adv+dCdtT_dis+r_T+kNAPL*(Csol_T-CsubT));
% Xylene
dCdt_X=(1./Rd_X).*(dCdtX_adv+dCdtX_dis+r_X+kNAPL*(Csol_X-CsubX));
% Mobile BIOMASS
dCdtBIO_m=dCdtBIO_adv+dCdtBIO_dis+rBIO_m-Cbio_m*b+Cbio_im*kdet*rho_bulk./poro-Cbio_m*katt;
% Immobile BIOMASS
dCdtBIO_im=rBIO_im-Cbio_im*b-Cbio_im*kdet+Cbio_m*katt.*poro/rho_bulk;
%Back to ODE form composed by both drawdown and concentration
dhdt=[dsdt;dCdtO;dCdt_B;dCdt_T;dCdt_X;dCdtBIO_m;dCdtBIO_im];
%Update waitbar
% waitbar(log10(t+1e-5)/log10(tMax),wait);

