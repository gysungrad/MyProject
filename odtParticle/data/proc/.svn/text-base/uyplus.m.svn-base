%
% Compute uplus versus yplus profile for channel flow cases
% Channel Flow Params: 
%
% Re_tau   = u_tau * delta / nu
% Re_tau   = delta/nu * sqrt(-delta/rho * dpdx)
% u_tau    = sqrt(-delta/rho * dpdx)
% delta_nu = nu/u_tau
% 
% delta is channel half width, nu is kinematic visc
%

clear;
clc;

%---------------------- DOLODT

utau     = 1.0; 
delta_nu = 2.53E-3/utau;     

data = load 'dolodt_Re395_retune';
data = load 'dolodt_Re395_retune2';

y = data(:,1);
u = data(:,3);

[npts di] = size(u);
i=1:npts/2;
j=npts:-1:npts/2+1;
y = y(1:npts/2);
u = (u(i) + u(j))/2;

up = u/utau           ;%    *sqrt(2); 
yp = y/delta_nu       ;%    *sqrt(2);

%---------------------- DNS

dnsdata = load 'dnsChannel_Re395means';
yp_dns1 = dnsdata(2:end,2);
up_dns1 = dnsdata(2:end,3);

dnsdata = load 'dnsChannel_Re590means';
yp_dns2 = dnsdata(2:end,2);
up_dns2 = dnsdata(2:end,3);

%---------------------- Bodt

bodt = load 'Bodt_Re395';
yp_bodt1 = bodt(2:end-1,1)/1.266213E-4;
up_bodt1 = bodt(2:end-1,2)/0.17692;

bodt = load 'Bodt_Re590';
yp_bodt2 = bodt(2:end-1,1)/8.4784E-5;
up_bodt2 = bodt(2:end-1,2)/0.17692;

[npts di] = size(up_bodt1);
i=1:npts/2;
j=npts:-1:npts/2+1;

yp_bodt1 =  yp_bodt1(1:npts/2);
up_bodt1 = (up_bodt1(i) + up_bodt1(j))/2;
yp_bodt2 =  yp_bodt2(1:npts/2);
up_bodt2 = (up_bodt2(i) + up_bodt2(j))/2;

%---------------------- ankodt

%ank = load 'ankodt';
ank = load 'ankodt_retune';
yp_nk1 = ank(:,1)/delta_nu ; %* sqrt(2);
up_nk1 = ank(:,2)/utau     ; %* sqrt(2); 

[npts di] = size(up_nk1);

i=1:npts/2;
j=npts:-1:npts/2+1;

yp_nk1 =  yp_nk1(1:npts/2);
up_nk1 = (up_nk1(i) + up_nk1(j))/2;

%---------------------- plot

semilogx(yp_dns1,up_dns1, '-', ...
         yp_bodt1,up_bodt1, '-', ...
%         yp_bodt2,up_bodt2, '-', ...
         yp_nk1,up_nk1, '-', ...
         yp,up,'-');
title('Channel Flow');
ylabel('u+');
xlabel('y+');
%legend('DNS 590', 'BODT 395', 'BODT 590', 'NKODT', 'DOL', 'location', 'west');
legend('DNS 395', 'BODT 395', 'NKODT', 'DOL', 'location', 'west');
axis([0.01 1000 0 25]);

