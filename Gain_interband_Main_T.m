%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 11February2020, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program computes the absorption and gain of interband transition in 
% semiconductor. It can computes the gain in bulk (3D) or in quantum well (2D)
% Various materail and their alloys are available in the library

% This special version allow to do a sweep over the Temperature T

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

T  = 250:10:350;                 % Temperature [K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N3d= 2e18*1e6;            % Carrier concentration [m-3] 
M  = InGaAs10;            % Choose the material from the library
%M  = GaAs;                % Choose the material from the library
d  = 2;                   % Dimension, bulk=3, Quantum well=2
Lqw= 10e-9;               % Quantum well width (meter) if d=2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EP  = M(idx_EP_K);        % EP Kane
me  = M(idx_me);          % electron mass
mhh = M(idx_mhh);         % heavy hole mass
mr  = me*mhh/(me+mhh);    % reduced mass
nopt=sqrt(M(idx_Epsi));   % optical index

FWHM=1e-2;                % homogeneous broadening (eV)

% in case we wanna have current instead of electron number
%I=10e-3;                 % injected current (A)
%S=pi*(10e-6)^2;          % surface pumped (m-2)
%tau=1e-9;                % total interband recombinason time (s)
%N3d=I/S/e/Lqw*tau        % injected charge [m-3]

N2d=N3d*Lqw;              % sheet density in case d=2 [m-2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=linspace(-0.25,2,1000);
EEx=repmat(E,length(E),1);
EEy=repmat(E',1,length(E));

%L=(FWHM/2)^2 * 1./ ( ( E-E' ).^2 + (FWHM/2).^2  ) ;
L=(FWHM/2)^2 * 1./ ( ( EEx-EEy ).^2 + (FWHM/2).^2  ) ;

%L=L./trapz(E,L,2);
L=L./repmat(trapz(E,L,2) , 1 ,length(E)  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(T)
  
Eg  = M(idx_Eg6c) - (M(idx_alphaG)*T(i)^2) ./ (T(i)+M(idx_betaG));   %Eg = Eg0 - (a*T.^2)./(T + b);
Echv = Eg + mr/me *(E-Eg);
Evhv =    - mr/mhh*(E-Eg);

if d==2
  [Efc,Efv,ron,rop,FEc,FEv,FEcc,FEvv,alpha,Gain]=Gain2D_interband_f(N2d,me,mhh,E,Eg,EP,L,T(i),nopt,FWHM,2,Lqw);
end
if d==3
  [Efc,Efv,ron,rop,FEc,FEv,FEcc,FEvv,alpha,Gain]=Gain3D_interband_f(N3d,me,mhh,E,Eg,EP,L,T(i),nopt,FWHM,3);
end

alphaT(:,i)=alpha';
GainT(:,i)=Gain;
s{i}=strcat('\fontsize{20}T=',num2str(T(i)),'K');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure('position',[-3500 10 1200 900])
figure('position',[100 100 1200 900])
FS=20;
c=colormap(jet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;grid on; box on;

for i=1:length(T)
  col=c((i-1)*round(64/length(T))+1,:);
  plot(E,GainT(:,i),'linewidth',2,'color',col)
end

xlabel('E (eV)')
ylabel('Losses / Gain (cm-1)')

xlim([0.9*Eg E(end)])
title(strcat('InGaAs-10%: ',num2str(d),'d; N=',num2str(N3d*1e-6),'cm-3'))
legend(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%