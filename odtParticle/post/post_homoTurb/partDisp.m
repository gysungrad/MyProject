%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guangyuan Sun 05/12
% Process homogeneous isotropic turbulence
% Standalone code.
% Mutiple realizations
% Calculate particle dispersion from particle dump files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function partDisp(dir, nrlzn, npart, jrf) 
% clc; clear;
% 
% dir = 'homoSL3_040714';
% nrlzn = 512;    % # of realization
% npart = 3;   % # of particles
ref = zeros(npart,1);
% jrf = 8; % # of reference 

%%%%%%%%%%%%%%%% read dumptime file %%%%%%%%%%%%%%%%%%%% 
command = ['../../input/dumpTimes.inp'];
file = fopen(command);
ln = fgetl(file);
clear dumpTime;
kk = 1;
while(~feof(file))
    ln = fgetl(file);
    dumpTime(kk,:) = [sscanf(ln, '%f')]';
    kk = kk+1;
end
fclose(file);

ndump = length(dumpTime);      % # of dumptimes
ref_time = dumpTime(jrf);      % reference time

dispersion = zeros(ndump,npart); % particle dispersion sqr root
dispersion_sqr = zeros(ndump,npart); %particle dipsersion sqr

for i = 1:nrlzn

%%%%%%%%%%%%%%%% calculate particle reference position %%%%%%%%%%%%%%%%%%%% 
    command = ['../../data/',dir, '/data_',num2str(i-1),'/dmp_part_',num2str(jrf),'.dat'];
    file = fopen(command);
    ln = fgetl(file);
    ln = fgetl(file);
    clear B;
    jj = 1;
    while(~feof(file))
        ln = fgetl(file);
        B(jj,:) = [sscanf(ln, '%f')]';
        jj = jj+1;
    end
    for kk = 1:npart
        ref(kk,1) = B(kk,2);
    end
    fclose(file);
    
%%%%%%%%%%%%%%%% calculate particle dipersion sqr %%%%%%%%%%%%%%%%%%%% 
    for j = 1:ndump
        command = ['../../data/', dir, '/data_',num2str(i-1),'/dmp_part_',num2str(j),'.dat'];
        file = fopen(command);
        ln = fgetl(file);
        ln = fgetl(file);
        clear A;
        ii = 1;
        while(~feof(file))
            ln = fgetl(file);
            A(ii,:) = [sscanf(ln, '%f')]';
            ii = ii+1;
        end
        dispersion_sqr(j,:) = dispersion_sqr(j,:) + ((A(:,2)-ref(:,1)).^2)'/nrlzn; % calculate dispersion sqr 
        fclose(file);
    end
end

dispersion = sqrt(dispersion_sqr);

%%%%%%%%%%%%%%% Open the output file for write %%%%%%%%%%%%%%%%%%%%%%%%%

outfilename=['partDisp_homogeneousSL_',dir,'.dat'];
[fiout,message]=fopen(outfilename,'w');
if (fiout<0)
    error(message);
end

%%%%%%%%%%%%%%% output the data %%%%%%%%%%%%%%%%%%%%%%%%%
for iDump=jrf:ndump
    fprintf(fiout,'%20.10E',dumpTime(iDump)-ref_time); % sec
    for iPart=1:npart
        fprintf(fiout,'%20.10E',dispersion_sqr(iDump,iPart)*10000); % cm^2
        fprintf(fiout,'%20.10E',dispersion(iDump,iPart)*100);       % cm
    end
    fprintf(fiout,'\n');   
end
fclose(fiout);

