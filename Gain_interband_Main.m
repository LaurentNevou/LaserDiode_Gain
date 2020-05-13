%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 11February2020, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program computes the absorption and gain of interband transition in 
% semiconductor. It can computes the gain in bulk (3D) or in quantum well (2D)
% Various materials and their alloys are available in the library

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library
TernaryAlloy;             % compute the ternary alloy
%QuaternaryAlloy;          % compute the quaternary alloy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N3d= 5e18*1e6;            % Carrier concentration [m-3] 
M  = InGaAs10;            % Choose the material from the library
%M  = GaAs;               % Choose the material from the library
T  = 300;                 % Temperature [K]
d  = 3;                   % Dimension, bulk=3, Quantum well=2
Lqw= 10e-9;               % Quantum well width (meter) if d=2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eg  = M(idx_Eg6c) - (M(idx_alphaG)*T^2) ./ (T+M(idx_betaG));   %Eg = Eg0 - (a*T.^2)./(T + b);
EP  = M(idx_EP_K);        % EP Kane
me  = M(idx_me);          % electron mass
mhh = M(idx_mhh);         % heavy hole mass
mr  = me*mhh/(me+mhh);    % reduced mass
nopt=sqrt(M(idx_Epsi));   % optical index

FWHM=1e-3;                % homogeneous broadening (eV)

% in case we wanna have current instead of electron number
%I=10e-3;                 % injected current (A)
%S=pi*(10e-6)^2;          % surface pumped (m-2)
%tau=1e-9;                % total interband recombinason time (s)
%N3d=I/S/e/Lqw*tau        % injected charge [m-3]

N2d=N3d*Lqw;              % sheet density in case d=2 [m-2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=linspace(-0.25,2,1000);
EEx=repmat(E,length(E),1);
EEy=repmat(E',1,length(E));

%L=(FWHM/2)^2 * 1./ ( ( E-E' ).^2 + (FWHM/2).^2  ) ;
L=(FWHM/2)^2 * 1./ ( ( EEx-EEy ).^2 + (FWHM/2).^2  ) ;

%L=L./trapz(E,L,2);
L=L./repmat(trapz(E,L,2) , 1 ,length(E)  );

Echv = Eg + mr/me *(E-Eg);
Evhv =    - mr/mhh*(E-Eg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if d==2
  [Efc,Efv,ron,rop,FEc,FEv,FEcc,FEvv,alpha,Gain]=Gain2D_interband_f(N2d,me,mhh,E,Eg,EP,L,T,nopt,FWHM,2,Lqw);
end
if d==3
  [Efc,Efv,ron,rop,FEc,FEv,FEcc,FEvv,alpha,Gain]=Gain3D_interband_f(N3d,me,mhh,E,Eg,EP,L,T,nopt,FWHM,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure('position',[-3500 10 1200 900])
figure('position',[100 100 1200 900])
FS=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1,'fontsize',FS)
hold on;grid on; box on;

plot(E,rop ,'b.-')
plot(E,rop.*FEv,'r.-')
plot(Efv*[1 1]   ,[0 1]*max(rop.*FEv)*5,'m')

xlim([-1 1]*0.25)
ylim([0 1]*max(rop.*FEv)*5)
xlabel('E (eV)')

if d==2
  ylabel('2d Density of state (eV-1.m-2)')
  title(strcat('2d: \color{red}Efv=', num2str(Efv*1e3,'%.1f'),'meV'))
elseif d==3
  ylabel('3d Density of state (eV-1.m-3)')
  title(strcat('3d: \color{red}Efv=', num2str(Efv*1e3,'%.1f'),'meV'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2,'fontsize',FS)
hold on;grid on; box on;

plot(E,ron ,'b.-')
plot(E,ron .* FEc,'g.-')
plot(Efc*[1 1]+Eg,[0 1]*max(ron.*FEc)*5,'m')

xlim([0.9 1.2]*Eg)
ylim([0 1]*max(ron.*FEc)*5)
xlabel('E (eV)')

if d==2
  ylabel('2d Density of state (eV-1.m-2)')
  title(strcat('2d: \color{green}Efc=',num2str(Efc*1e3,'%.1f'),'meV'))
elseif d==3
  ylabel('3d Density of state (eV-1.m-3)')
  title(strcat('3d: \color{green}Efc=',num2str(Efc*1e3,'%.1f'),'meV'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3,'fontsize',FS)
hold on; grid on;box on;

plot(E,FEcc,'g.-')
plot(E,FEvv,'r.-')
plot(E,FEcc-FEvv,'b.-')

xlabel('E (eV)')
ylabel('Losses / Gain')

xlim([0 E(end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4,'fontsize',FS)
hold on; grid on;box on;

plot(E,Gain,'m.-')

xlabel('E (eV)')
ylabel('Losses / Gain (cm-1)')

xlim([0.9*Eg E(end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
