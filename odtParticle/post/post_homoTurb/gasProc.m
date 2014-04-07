

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process homogeneous isotropic turbulence
% Standalone code.
% Mutiple realizations
% Process gas phase data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gasProc(dir, Nrlzn, a, b, U)
%clc; clear;

%dir = 'homoSL1_040614';
% dir = 'SLHollowI_1023';
% Nrlzn = 256;   % number of realizations
nx = 150;      % number of interpolation points
init_pt = 75;  % number of the point which we calculate corr coefficients from
% a = 0;         % left endpoint of the domain
% b = .5;        % right endpoint of the domain
visc = 1.8E-5/1.2; % viscosity
iuvel = 6;     % index of u velocity in dump files


M = .0254; % Size of mesh
% U = 6.55; % Bulk streamwise velocity

%%%%%%%%%%%%%%% get dump times %%%%%%%%%%%%%%%%

command = ['../../input/dumpTimes.inp']
file = fopen(command);
ln = fgetl(file);
k = 1;
while(~feof(file))
    ln = fgetl(file);
    dataTime(k,:) = [sscanf(ln, '%f')]';
    k = k+1;
end
%ndump = length(dataTime) - 1;
ndump = length(dataTime);
fclose(file);    

% Get streamwise positions at dumptimes and scale by mesh length
dump_pos = dataTime*U;
scaled_dump_pos = dump_pos/M;

% Set up vector for interpolation

X = linspace(a,b,nx)';
dx = (b-a)/nx;

% Initialize relevant variables

varianceU = zeros(ndump,1);
autocorr = zeros(nx,ndump);
spectrum = zeros(nx,ndump);

dissip = zeros(ndump,1);
% Loop over all dump times
for j = 1:ndump
    
    ensemble_mean_U_sqr = zeros(nx,1);
    % Loop over all realizations
    for ifile = 1:Nrlzn
        [ifile j]
        % Read in data from dump file
        command = ['../../data/',dir,'/data_',num2str(ifile-1),'/dmp_odtl_',num2str(j),'.dat']
        file = fopen(command);

        ln = fgetl(file);
        ln = fgetl(file);
        ln = fgetl(file);
        ln = fgetl(file);

        clear A;
        i = 1;
        while(~feof(file))
            ln = fgetl(file);
            A(i,:) = [sscanf(ln, '%f')]';
            i = i+1;
        end
        
        % Interpolate u velocity on the line
        interp_u = interp1(A(:,1),A(:,6),X,'spline','extrap'); 
        
        % Calculate the u velocity derivative in the line direction
        N = length(A(:,1));
        deriv = zeros(N,1);
        deriv(1) = (A(2,6)-A(1,6))/(A(2,1)-A(1,1));
        deriv(N) = (A(N,6)-A(N-1,6))/(A(N,1)-A(N-1,1));
        for i = 2:N-1
            deriv(i) = (A(i+1,6)-A(i,6))/(A(i+1,1)-A(i,1));
        end
        
        % Interpolate the derivative
        interp_der = interp1(A(:,1),deriv,X,'spline','extrap');

        mean_sqrU = mean(interp_u.^2);  % mean of velocity square
        mean_U = mean(interp_u);     % mean of velocity
        sqr_mean_U = mean_U^2; % square of mean velocity
        varianceU(j,1) = varianceU(j,1) + (mean_sqrU - sqr_mean_U)/Nrlzn;
        
        % Add up terms for numerator of autocorrelation
        autocorr(:,j) = autocorr(:,j) + interp_u(init_pt)*interp_u/Nrlzn;
        % Add up terms for the denominator of autocorrelation
        ensemble_mean_U_sqr = ensemble_mean_U_sqr + interp_u.^2/Nrlzn;
        
        % Add up derivative terms for dissipation rate
        dissip(j) = dissip(j) + mean(interp_der.^2)/Nrlzn;
        
        fclose(file);
    end
    
    % Divide for autocorrelation
    autocorr(:,j) = autocorr(:,j)/ensemble_mean_U_sqr(init_pt);
    % Take fft to calculate the spectrum and scale it
    spectrum(:,j) = fft(autocorr(:,j))/nx;
end

% Calculate dissipation rate and kinetic energy
dissip = 3*visc*dissip;
kin_ener = 3/2*varianceU;

% Compute and scale like Snyder Lumley or Wells Stock
%decay = U^2./(varianceU.^2)/1000;
decay = sqrt(varianceU)/U*100;
norm_mean = sqrt(2/3*kin_ener)/U*100;
integral_time = kin_ener./dissip*1000;
integral_length = kin_ener.^(3/2)./dissip*100;
kolm_time = sqrt(visc./dissip)*1000;
kolm_length = (visc^3./dissip).^(1/4)*100;
dissip = dissip*10000;

% Output the data

% var = [dataTime(2:ndump+1) dissip kin_ener decay integral_time integral_length kolm_time kolm_length];
var = [scaled_dump_pos(1:ndump) dissip kin_ener decay integral_time integral_length kolm_time kolm_length];
filename = ['GasSL',dir,'.dat'];
%filename = ['Gas_Time_Scales_',dir,'.dat'];
save(filename, 'var', '-ascii');

var = [X autocorr];
filename = ['Gas_auto_',dir,'.dat'];
save(filename, 'var', '-ascii');

k = (0:nx-1)'/(b-a);
var = [k abs(spectrum)];
filename = ['Gas_spect_',dir,'.dat'];
save(filename, 'var', '-ascii');



