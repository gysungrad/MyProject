% run as standalone (uncomment as noted), or use pdfdrv.m driver
%
% Computes the PDF of accepted eddies
% Read in a file called out that has sizes of each accepted eddy as the columnPos column.
% Change nbins, and columnPos
% Plots the pdf with the model pdf and the model pdf with a new beta
%


%clc;
%clear;

nbins = 60;
columnPos = 6;

%---------------- compute eddy pdf

%oo = load 'out_re395_1';                       % uncomment standalone

dat = oo(:,columnPos);
[npts id] = size(dat);

Lmin = 0.05*min(dat);
Lmax = 2.3*max(dat);
Lmin = 0.001;
Lmax = 10;

binf = logspace(log10(Lmin), log10(Lmax), nbins+1);
pdf = zeros(nbins,1);

for i=1:nbins

    jpos = find( dat >= binf(i) & dat < binf(i+1) );
    [ni id] = size(jpos);
    pdf(i) = pdf(i) + ni;

end

pdf = pdf / npts;
ds = zeros(nbins,1);
for i=1:nbins
    ds(i) = binf(i+1) - binf(i);
    pdf(i) = pdf(i) / ds(i);
end
ds = ds';

i=1:nbins;
bins = (binf(i+1) + binf(i))/2;

%--------------- compare to drawn distribution

%domainLength = 2;                              % uncomment standalone
Lp   = 0.015 * domainLength;
Lmax = 1.0   * domainLength;
Lmin = 0.004 * domainLength;
Cmax = exp(-2*Lp/Lmax);
Cmin = exp(-2*Lp/Lmin);
C    = 2*Lp/(Cmax-Cmin);
f    = C./bins.^2 .*exp(-2*Lp./bins);
f = f';

%beta  = -2.5;
%gamma = -beta/(beta+1);
%C     = 0.0429;              % trial and error to get ss=1 
%f2    = C.*bins.^beta .* exp(gamma*(Lp./bins).^-(beta+1));
%ss = sum(f2'.*ds)
%Lp = bins(find(pdf==max(pdf)))

%--------------- draw from the distribution to verify the results

%nSe = 60000;         % number of eddies to draw
%rr = rand(nSe,1);
%Seddy = -2.0*Lp./log(rr*(Cmax-Cmin)+Cmin);  % eddy distribution
%pdfSe = zeros(nbins,1);                     % pdf of eddies
%ntot = 0;
%%for i=1:nbins
%    jpos = find( Seddy >= binf(i) & Seddy < binf(i+1) );
%    [ni id] = size(jpos);
%    pdfSe(i) = pdfSe(i) + ni;
%end
%pdfSe = pdfSe/nSe;
%for i=1:nbins
%    pdfSe(i) = pdfSe(i) / ds(i);
%end


%--------------- plot the results
   
%semilogx((bins),pdf,'o-', bins,f,'-', bins,f2,'-');
%semilogx((bins),pdf,'o-', bins,f,'-', bins,pdfSe);
semilogx((bins),pdf,'o-', bins,f);
xlabel('Eddy Size');
ylabel('PDF');
title('PDF of Accepted Eddies');
%legend('Eddy PDF', 'Model PDF', 'Model PDF (new beta)');
%legend('Eddy PDF', 'Model PDF', 'sample MPDF');
legend('Eddy PDF', 'Model PDF'); 
axis([0.5*bins(1) 1.5*bins(end) f(1) 1.2*max([pdf ;f])]);


