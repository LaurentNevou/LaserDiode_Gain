function[Efc_2d,Efv_2d,ro2dn,ro2dp,FEc_2d,FEv_2d,FEcc_2d,FEvv_2d,alpha2D,Gain2D]=Gain2D_interband_f(N2d,me,mhh,E,Eg,EP,L,T,nopt,FWHM,d,Lqw)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rosencher, "Optoelectronics", 
% chap7: "Optical properties of semiconductor"
% 7.3 Optical susceptibility of a semiconductor
% 7.6 Conditions for optical amplification in semiconductors
% chap13: "Light emitting diodes and laser diodes"
% 13.4 Optical amplification in heterojunction diodes 
% 13.6 Quantum well laser diode

% J. Faist, "Optical properties of semiconductor"
% chap5: "Bulk semiconductors: bandstructure and fundamental gap"
% 5.2 Computation of the absorption edge in bulk materails
% 5.4 Gain

% S. L. Chuang: "Physics os Optoelectronic Devices"
% Part 3: Generation of light
% chap9: Optical Processes in semiconductor

% G. Bastard "Wave mechanics applied to semiconductor heterostructures"
% chap7 "Optical properties of quasi bi-dimensional systems", page 241

% G. Fishman "Semi-conducteurs: les bases de la theorie k.p"
% https://www.eyrolles.com/Sciences/Livre/semi-conducteurs-les-bases-de-la-theorie-k-p-9782730214971/
% https://www.editions.polytechnique.fr/?afficherfiche=149
% chap10: Absorption
% 10.2.2 Coefficiant d'absorption

% https://www.nextnano.de/nextnano3/tutorial/1Dtutorial12.htm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h    = 6.62606896E-34;               %% Planck constant J.s
hbar = h/(2*pi);
e    = 1.602176487E-19;              %% charge de l electron Coulomb
c    = 2.99792458E8;                 %% vitesse de la lumiere m/s
Epsi0= 8.854187817620E-12;           %% constant dielectric du vide F/m
kB   = 1.3806488E-23;                %% Boltzmann's constant (J/K)
m0   = 9.10938188E-31;               %% electron mass kg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mr=me*mhh/(me+mhh);
Echv = Eg + mr/me *(E-Eg);
Evhv =    - mr/mhh*(E-Eg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Efc_2d = +Fermi2D_f(N2d,me ,T,FWHM);
Efv_2d = -Fermi2D_f(N2d,mhh,T,FWHM);

ro2dn  = 0*E;
ro2dp  = 0*E;
ro2dEg = 0*E;

ro2dn(E>Eg)  = e*me*m0 /(pi*(hbar)^2);
ro2dp(E<0)   = e*mhh*m0/(pi*(hbar)^2);
ro2dEg(E>Eg) = e*mr*m0 /(pi*(hbar)^2);
ro2dn  = trapz(E, repmat(ro2dn ,[length(E) 1]) .* L , 2)';
ro2dp  = trapz(E, repmat(ro2dp ,[length(E) 1]) .* L , 2)';
ro2dEg = trapz(E, repmat(ro2dEg,[length(E) 1]) .* L , 2)';

FEc_2d = 1./(1+exp( ( E-Efc_2d -Eg )/(kB*T/e))) ; 
FEv_2d = 1./(1+exp(-( E-Efv_2d     )/(kB*T/e))) ;

FEcc_2d = 1./(1+exp(( Echv-Efc_2d -Eg )/(kB*T/e))) ;
FEvv_2d = 1./(1+exp(( Evhv-Efv_2d     )/(kB*T/e))) ; 

Ne_2d = trapz(E,ro2dn .* FEc_2d);
Nh_2d = trapz(E,ro2dp .* FEv_2d);
N2d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha2D=2*pi^2*e^2/(6*pi*Epsi0*nopt*m0*c) * (EP/Eg) * ro2dEg ./ (e*Eg/hbar) / Lqw; % [m-1]
Gain2D = (alpha2D*1e-2.*(FEcc_2d-FEvv_2d))'; % [cm-1]