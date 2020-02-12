%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% III-V Ternary alloys on InP %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.53;
InGaAs = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;   % In0.53Ga0.47As is lattice matched on InP
x=0.37;
InGaAs37 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;   % In0.37Ga0.63As STRAINED on InP
x=0.45;
InGaAs45 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;   % In0.45Ga0.55As STRAINED on InP
x=0.715;
InGaAs71 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;   % In0.715Ga0.285As STRAINED on InP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.48;
AlInAs = x*AlAs + (1-x)*InAs - x*(1-x)*AlInAs_bowing;   % Al0.48In0.52As is lattice matched on InP
x=0.80;
AlInAs80 = x*AlAs + (1-x)*InAs - x*(1-x)*AlInAs_bowing;   % Al0.80In0.20As STRAINED on InP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.51;
GaAsSb = x*GaAs + (1-x)*GaSb - x*(1-x)*GaAsSb_bowing;   % GaAs0.51Sb0.49 is lattice matched on InP
% QCL: Al-fee by G. Strasser's group in Vienna, JVSTB, 32. 02C104 (2014)
% I did adjust VBO_bowing in order to get CBO=0.36eV between InGaAs/GaAsSb

x=0.56;
AlAsSb = x*AlAs + (1-x)*AlSb - x*(1-x)*AlAsSb_bowing;   % AlAs0.56Sb0.44 is lattice matched on InP
% the X-valley is 0.5eV above Gamma point (indirect bandgap material)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% III-V Ternary alloys on GaAs %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.10;
AlGaAs_bowing(idx_Eg6c) = -0.127 + 1.310*x;   % from Vurgaftman
AlGaAs_bowing(idx_VBO)  =  0.27  -  1.05*x;   % to adjust VBO in order to fit CBO from THz QCL
AlGaAs10 = x*AlAs + (1-x)*GaAs - x*(1-x)*AlGaAs_bowing;

x=0.15;
AlGaAs_bowing(idx_Eg6c) = -0.127 + 1.310*x;   % from Vurgaftman
AlGaAs_bowing(idx_VBO)  =  0.27  -  1.05*x;   % to adjust VBO in order to fit CBO from THz QCL
AlGaAs15 = x*AlAs + (1-x)*GaAs - x*(1-x)*AlGaAs_bowing;

x=0.20;
AlGaAs_bowing(idx_Eg6c) = -0.127 + 1.310*x;   % from Vurgaftman
AlGaAs_bowing(idx_VBO)  =  0.27  -  1.05*x;   % to adjust VBO in order to fit CBO from THz QCL
AlGaAs20 = x*AlAs + (1-x)*GaAs - x*(1-x)*AlGaAs_bowing;

x=0.30;
AlGaAs_bowing(idx_Eg6c) = -0.127 + 1.310*x;   % from Vurgaftman
AlGaAs_bowing(idx_VBO)  =  0.27  -  1.05*x;   % to adjust VBO in order to fit CBO from THz QCL
AlGaAs30 = x*AlAs + (1-x)*GaAs - x*(1-x)*AlGaAs_bowing;

x=0.40;
AlGaAs_bowing(idx_Eg6c) = -0.127 + 1.310*x;   % from Vurgaftman
AlGaAs_bowing(idx_VBO)  =  0.27  -  1.05*x;   % to adjust VBO in order to fit CBO from THz QCL
AlGaAs40 = x*AlAs + (1-x)*GaAs - x*(1-x)*AlGaAs_bowing;

x=0.50;
AlGaAs_bowing(idx_Eg6c) = -0.127 + 1.310*x;   % from Vurgaftman
AlGaAs_bowing(idx_VBO)  =  0.27  -  1.05*x;   % to adjust VBO in order to fit CBO from THz QCL
AlGaAs50 = x*AlAs + (1-x)*GaAs - x*(1-x)*AlGaAs_bowing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.05;
InGaAs05 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.06;
InGaAs06 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.07;
InGaAs07 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.08;
InGaAs08 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.09;
InGaAs09 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.10;
InGaAs10 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.15;
InGaAs15 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.20;
InGaAs20 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.25;
InGaAs25 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.30;
InGaAs30 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;
x=0.40;
InGaAs40 = x*InAs + (1-x)*GaAs - x*(1-x)*InGaAs_bowing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.49;
GaInP = x*InP + (1-x)*GaP - x*(1-x)*GaInP_bowing;   % Ga0.51In0.49P is lattice matched on GaAs

x=0.48;
AlInP = x*InP + (1-x)*AlP - x*(1-x)*AlInP_bowing;   % Al0.52In0.48P is lattice matched on GaAs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% III-V Ternary alloys on InAs/GaSb/AlSb %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.25;
GaInSb25 =  x*InSb + (1-x)*GaSb - x*(1-x)*GaInSb_bowing;
x=0.35;
GaInSb35 =  x*InSb + (1-x)*GaSb - x*(1-x)*GaInSb_bowing;

x=0.09;
GaAsSb91 = x*GaAs + (1-x)*GaSb - x*(1-x)*GaAsSb_bowing;   % GaAs0.09Sb0.91 is lattice matched on InAs

x=0.91;
InAsSb   = x*InAs + (1-x)*InSb - x.*(1-x)*InAsSb_bowing;    % In0.09As0.91Sb is lattice matched on GaSb, Eg~0.3eV

x=0.5;y=0.08;
AlGaAsSb = x*GaSb + (1-x)*(y*AlAs + (1-y)*AlSb - y*(1-y)*AlAsSb_bowing); % high barrier lattice matched on GaSb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% IV-IV Ternary alloys %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0.25;
SiGe25 = x*Ge + (1-x)*Si - x*(1-x)*SiGe_bowing;
x=0.3;
SiGe30 = x*Ge + (1-x)*Si - x*(1-x)*SiGe_bowing;
x=0.35;
SiGe35 = x*Ge + (1-x)*Si - x*(1-x)*SiGe_bowing;
x=0.4;
SiGe40 = x*Ge + (1-x)*Si - x*(1-x)*SiGe_bowing;
x=0.5;
SiGe50 = x*Ge + (1-x)*Si - x*(1-x)*SiGe_bowing;
x=0.8;
SiGe80 = x*Ge + (1-x)*Si - x*(1-x)*SiGe_bowing;
