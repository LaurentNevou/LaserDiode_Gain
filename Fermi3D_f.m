function[Ef]=Fermi3D_f(N3d,meff,T,FWHM)

de=0.001;
E1=linspace( -2         , -0.1 , 500  );
E2=linspace( E1(end)+de ,  0   , 500  );
E3=linspace( E2(end)+de ,  0.5 , 5000 );
E4=linspace( E3(end)+de ,  5   , 500  );

E=sort([E1 E2 E3 E4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant J.s
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron Coulomb
kB =1.3806488E-23;              %% Boltzmann's constant (J/K)
m0=9.10938188E-31;              %% electron mass kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% density of state %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%L=(FWHM/2)^2 * 1./ ( ( E-E' ).^2 + (FWHM/2).^2  ) ;
%L=L./trapz(E,L,2);

ro3d=0*E;
cst=1/(2*pi^2) *  (2*e*meff *m0/(hbar^2))^(3/2);
ro3d(E>0) = cst * sqrt(E(E>0));  
%ro3d=trapz(E, repmat(ro3d,[length(E) 1]) .* L , 2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% calcul of Fermi level at T=0K %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ef0=(3/2*N3d/cst)^(2/3);  % This is not very correct since I had a broadening in ro(E)...

Ef=Ef0;
FE = 1./(1+exp((E-Ef)/(kB*T/e))) ;
N3dX= trapz(E,ro3d.*FE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% calcul of Fermi level at any temperature %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dE=0.025;   % step of the scan in eV to pass below the Fermi level
s=1;

while  abs(N3dX - N3d)/N3d > 1e-10

  if  N3dX > N3d
    if s==0
      dE=dE/2;
    end
    s=1;
    Ef= Ef - dE;    
    FE= 1./(1+exp((E-Ef)/(kB*T/e))) ; 
    N3dX = trapz(E,ro3d.*FE);
    
  else
    if s==1
      dE=dE/2;
    end
    s=0;
    Ef= Ef + dE;    
    FE= 1./(1+exp((E-Ef)/(kB*T/e))) ; 
    N3dX = trapz(E,ro3d.*FE);
  
  end
    
end       



%end