%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Extract general parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lattice parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'a');
%  if idx==1
%    idx_a=i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epsilon static

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Epsi_stat');
  if idx==1
    idx_Epsi=i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP from the Kane model

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'EP_Kane');
  if idx==1
    idx_EP_K = i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP from the Luttinger model

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'EP_Luttinger');
%  if idx==1
%    idx_EP_L = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 7c

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'Eg7c');
%  if idx==1
%    idx_Eg7c = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 8c

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'Eg8c');
%  if idx==1
%    idx_Eg8c = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 6v

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'Eg6v');
%  if idx==1
%    idx_Eg6v = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 6c (the band gap in a direct band gap semiconductor)

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Eg6cG');
  if idx==1
    idx_Eg6c = i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "alpha" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'alphaG');
  if idx==1
    idx_alphaG = i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "beta" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'betaG');
  if idx==1
    idx_betaG = i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split-off band energy Delta-so

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'Dso');
%  if idx==1
%    idx_Dso = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter for the electrons

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'F');
%  if idx==1
%    idx_F = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma1 

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'gamma1');
%  if idx==1
%    idx_g1 = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma2

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'gamma2');
%  if idx==1
%    idx_g2 = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma3

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'gamma3');
%  if idx==1
%    idx_g3 = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Valence band offset

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'VBO');
  if idx==1
    idx_VBO = i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mass for the 1band Schrodinger solver

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'me');
  if idx==1
    idx_me = i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mass for the 1band Schrodinger solver

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'mhh');
  if idx==1
    idx_mhh = i-1;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Strain parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ac Conduction band Strain parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'ac');
%  if idx==1
%    idx_ac = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% av Valence band Strain parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'av');
%  if idx==1
%    idx_av = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bv Valence band Strain parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'bv');
%  if idx==1
%    idx_bv = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dv Valence band Strain parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'dv');
%  if idx==1
%    idx_dv = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c11 Valence band Strain parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'c11');
%  if idx==1
%    idx_c11 = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c12 Valence band Strain parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'c12');
%  if idx==1
%    idx_c12 = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c44 Valence band Strain parameter

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'c44');
%  if idx==1
%    idx_c44 = i-1;
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%