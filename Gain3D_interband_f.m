function[Efc_3d,Efv_3d,ro3dn,ro3dp,FEc_3d,FEv_3d,FEcc_3d,FEvv_3d,alpha3D,Gain3D]=Gain3D_interband_f(N3d,me,mhh,E,Eg,EP,L,T,nopt,FWHM,d)

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

Efc_3d = +Fermi3D_f(N3d,me ,T,FWHM);
Efv_3d = -Fermi3D_f(N3d,mhh,T,FWHM);

ro3dn  = 0*E;
ro3dp  = 0*E;
ro3dEg = 0*E;

ro3dn(E>Eg)  = 1/(2*pi^2) *  (2*e*me *m0/(hbar^2))^(3/2) * sqrt(E(E>Eg)-Eg);
ro3dp(E<0  ) = 1/(2*pi^2) *  (2*e*mhh*m0/(hbar^2))^(3/2) * sqrt(-E(E<0));
ro3dEg(E>Eg) = 1/(2*pi^2) *  (2*e*mr *m0/(hbar^2))^(3/2) * sqrt(E(E>Eg)-Eg);
ro3dn  = trapz(E, repmat(ro3dn ,[length(E) 1]) .* L , 2)';
ro3dp  = trapz(E, repmat(ro3dp ,[length(E) 1]) .* L , 2)';
ro3dEg = trapz(E, repmat(ro3dEg,[length(E) 1]) .* L , 2)';

FEc_3d = 1./(1+exp( ( E-Efc_3d -Eg )/(kB*T/e))) ; 
FEv_3d = 1./(1+exp(-( E-Efv_3d     )/(kB*T/e))) ;

FEcc_3d = 1./(1+exp(( Echv-Efc_3d -Eg )/(kB*T/e))) ;
FEvv_3d = 1./(1+exp(( Evhv-Efv_3d     )/(kB*T/e))) ; 

Ne_3d = trapz(E,ro3dn .* FEc_3d);
Nh_3d = trapz(E,ro3dp .* FEv_3d);
N3d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha3D=2*pi^2*e^2/(6*pi*Epsi0*nopt*m0*c) * (EP/Eg) * ro3dEg ./ (e*Eg/hbar); % [m-1]
Gain3D = (alpha3D*1e-2.*(FEcc_3d-FEvv_3d))'; % [cm-1]