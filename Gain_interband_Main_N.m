%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 11February2020, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program computes the absorption and gain of interband transition in 
% semiconductor. It can computes the gain in bulk (3D) or in quantum well (2D)
% Various material and their alloys are available in the library

% This special version allow to do a sweep over the density of charges N

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

N3d= [1:5]*1e18*1e6;      % Carrier concentration [m-3] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M  = InGaAs10;            % Choose the material from the library
%M  = GaAs;                % Choose the material from the library
T  = 300;                 % Temperature [K]
d  = 2;                   % Dimension, bulk=3, Quantum well=2
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

Echv = Eg + mr/me *(E-Eg);
Evhv =    - mr/mhh*(E-Eg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(N3d)

if d==2
  [Efc,Efv,ron,rop,FEc,FEv,FEcc,FEvv,alpha,Gain]=Gain2D_interband_f(N2d(i),me,mhh,E,Eg,EP,L,T,nopt,FWHM,2,Lqw);
end
if d==3
  [Efc,Efv,ron,rop,FEc,FEv,FEcc,FEvv,alpha,Gain]=Gain3D_interband_f(N3d(i),me,mhh,E,Eg,EP,L,T,nopt,FWHM,3);
end

alphaN(:,i)=alpha';
GainN(:,i)=Gain;
s{i}=strcat('\fontsize{20}N=',num2str(N3d(i)*1e-6),'cm-3');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;grid on; box on;


plot(E,GainN,'linewidth',2)

xlabel('E (eV)')
ylabel('Losses / Gain (cm-1)')

xlim([0.9*Eg E(end)])
title(strcat(num2str(d),'d; T=',num2str(T),'K'))
legend(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
