%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Library load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DB = importdata('materialDB_ZB.csv',',');
%DB = importdata('materialDB_WZ.csv',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Here is a patch to be able to load the tables in Matlab AND Octave %%%%%%
% Matlab see the header in multiple cells while Octave see the header in one cell only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaAs';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaAs = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaAs=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlAs';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlAs = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlAs=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='InAs';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    InAs = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        InAs=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='InSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    InSb = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        InSb=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaSb = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaSb=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlSb = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlSb=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='InP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    InP = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        InP=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaP = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaP=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlP = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlP=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='Si';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    Si = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        Si=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='Ge';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    Ge = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        Ge=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='SiGe';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    SiGe_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        SiGe_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlGaAs';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlGaAs_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlGaAs_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='InGaAs';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    InGaAs_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        InGaAs_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlInAs';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlInAs_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlInAs_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaInP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaInP_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaInP_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlInP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlInP_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlInP_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlGaP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlGaP_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlGaP_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaInSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaInSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaInSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlInSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlInSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlInSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlGaSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlGaSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlGaSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaAsSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaAsSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaAsSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='InAsSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    InAsSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        InAsSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlAsSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlAsSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlAsSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaAsP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaAsP_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaAsP_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='InAsP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    InAsP_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        InAsP_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlAsP';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlAsP_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlAsP_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='GaPSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    GaPSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        GaPSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='InPSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    InPSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        InPSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='AlPSb';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    AlPSb_bowing = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        AlPSb_bowing=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%