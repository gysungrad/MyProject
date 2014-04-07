
clc;
clear;

load 'eddy_init.dat';
X1a = eddy_init;

u = X1a(:,5);
p = X1a(:,1);
pf = X1a(:,2);
pEnd = p(end) + p(end)-pf(end);

n = size(p);
n = n(1,1);
nf = n+1;
n2 = 2*n;

u2 = zeros(n2,1);
i = 1:2:n2;
u2(i)   = u;
u2(i+1) = u;

pf2 = zeros(n2,1);
pf2(i) = pf;
pf2(i+1) = [pf(2:end);pEnd];


plot(pf2, u2, '+-', p, u, 'o');
%axis([0 0.2 min(u)-0.1*abs(min(u)) max(u)+0.1*abs(max(u))]);
%axis([1.8 2 min(u)-0.1*abs(min(u)) max(u)+0.1*abs(max(u))]);

