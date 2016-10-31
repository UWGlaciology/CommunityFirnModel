
rhoi=917.0;
rhoig=rhoi/1000;
rho=350:rhoi;
rhog=rho/1000;
T1=265;
T2=220;
k_ice1=9.828*exp(-5.7e-3*T1);
k_ice2=9.828*exp(-5.7e-3*T2);

k1a=0.0206 + 0.7828*rhog.^2 + 2.472*rhog.^4; %Van Dusen (1929) from Cuffey and Paterson

k1b=2*k_ice1*rhog./(3*rhoig-rhog); % Schwerdfeger (1963), from C&P
k1c=2*k_ice2*rhog./(3*rhoig-rhog); % Schwerdfeger (1963), from C&P

k2=0.029*(1+100*rhog.^2); %Fujita, 2000

k3=k_ice*(rho./rhoi).^(2-0.5*rho/rhoi); %Schwander, 1997

k4=0.138-1.01*rhog+3.233*rhog.^2; %Sturm, 1997 (snow)

fig3=figure(3);
clf;
hold on;
box on;
grid on;
plot(rho,k1a,'b','linewidth',1.5)
plot(rho,k1b,'r','linewidth',1.5)
plot(rho,k1c,'r--','linewidth',1.5)
plot(rho,k2,'k','linewidth',1.5)
plot(rho,k3,'g','linewidth',1.5)
plot(rho,k4,'m','linewidth',1.5)
set(gca,'fontsize',14)
xlabel('Density (kg m^{-3})','fontsize',14)
ylabel('Thermal Conductivity (W  m^{-1} K^{-1})','fontsize',14)
leg1=legend('Van Dusen (1929)','Schwerdfeger (1963), 265K','Schwerdfeger (1963), 220K','Fujita (2000)','Schwander (1997)','Sturm (1997)');
set(leg1,'location','best','fontsize',12)
xlim([350 917])
savename='firn_thermal.eps';
print(fig3,savename,'-depsc');
