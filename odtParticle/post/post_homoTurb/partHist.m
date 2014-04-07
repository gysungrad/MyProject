%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guangyuan Sun 05/12
% Process homogeneous isotropic turbulence
% Standalone code.
% Mutiple realizations
% Calculate particle line velocity and dispersion from particle history file 
% Note that each particle history file in different rlzn may have different dumptime steps
%           So do the interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function partHist(dir, nrlzn, npart, timeDomain, ntime, Uinit)
% clc; clear;
% 
% dir = 'homoSL3_040714';
% nrlzn = 512;    % # of realization
% npart = 3;      % # of particles
% Uinit = 6.55;   % initial gas u velocity 
% ntime = 10000;  % # of interpolation time
% timeDomain = 1.6; % running time
time = linspace(0, timeDomain, ntime);

mean_vvel = zeros(ntime,npart); % particle line velocity
fluct_vvel_sqr = zeros(ntime,npart); % particle line velocity fluctuation sqr
fluct_vvel = zeros(ntime,npart); % particle line velocity fluctuation
dispersion = zeros(ntime,npart); % particle dispersion sqr root
dispersion_sqr = zeros(ntime,npart); %particle dipsersion sqr

for i = 1:nrlzn

    command = ['../../data/', dir, '/data_',num2str(i-1),'/hist_part.dat']
    file = fopen(command);
    ln = fgetl(file);

    clear A;
    ii = 1;
    while(~feof(file))
        ln = fgetl(file);
        A(ii,:) = [sscanf(ln, '%f')]';
        ii = ii+1;
    end

%%%%%%%%%%% interpolate particle vvel and dispersion in history files %%%%%%%%%%%%%%%%

    histTime = A(:,1);
    [histTime, index] = sort(histTime);
    uniq = [true, diff(histTime') ~= 0];
    for iPart = 1:npart
        histVvel = A(:, iPart*2+1);
        histVvel = histVvel(index);
        vvel = interp1(histTime(uniq), histVvel(uniq), time', 'spline'); 
        histDisp = A(:, iPart*2);
        histDisp = histDisp(index);
        disp = interp1(histTime(uniq), histDisp(uniq), time', 'spline');

        mean_vvel(:,iPart) = mean_vvel(:,iPart) + vvel/nrlzn;
        fluct_vvel_sqr(:,iPart) = fluct_vvel_sqr(:,iPart) + vvel.^2 / nrlzn; 
        dispersion_sqr(:,iPart) = dispersion_sqr(:,iPart) + disp.^2 / nrlzn; % calculate dispersion sqr 
    end
    fclose(file);
end


fluct_vvel_sqr = fluct_vvel_sqr - mean_vvel.^2;
fluct_vvel = sqrt(fluct_vvel_sqr);
dispersion = sqrt(dispersion_sqr);


%%%%%%%%%%%%%%% Open the output file for write %%%%%%%%%%%%%%%%%%%%%%%%%

outfilename=['partHist_homogeneousSL_',dir,'.dat'];
[fiout,message]=fopen(outfilename,'w');
if (fiout<0)
    error(message);
end

%%%%%%%%%%%%%%% output the data %%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:ntime
    fprintf(fiout,'%20.10E',time(m)-0.2652);
    for ipart=1:npart
        fprintf(fiout,'%20.10E',Uinit*Uinit/fluct_vvel_sqr(m,ipart));
        fprintf(fiout,'%20.10E',dispersion_sqr(m,ipart)*10000); % cm^2
    end
    fprintf(fiout,'\n');   
end

fclose(fiout);

