% Octave (Matlab) code to output cell sizes and velocity differences

clc;
clear;

data = load '10z.dat';

pf = data(:,2);
u  = data(:,5);

n = size(pf); n=n(1,1);
i = 1:n-1;
dx = pf(i+1)-pf(i);
du = u(i+1) - u(i);

subplot(3,1,1);
plot(pf,u,'x-');

subplot(3,1,2);
plot(pf(1:end-1),dx);

subplot(3,1,3);
plot(pf(1:end-1),du);


