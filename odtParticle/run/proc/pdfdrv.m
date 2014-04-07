% For plotting the accepted eddy pdf

clc;
clear;

domainLength = 2;
oo = load 'out_re395_1';
Pdf;
p1 = pdf;
f1 = f;

domainLength = 1;
oo = load 'out_re395_2';
Pdf;
p2 = pdf;
f2 = f;

domainLength = 0.5;
oo = load 'out_re395_3';
Pdf;
p3 = pdf;
f3 = f;

domainLength = 2;
oo = load 'out_re395_4';
Pdf;
p4 = pdf;
f4 = f;

domainLength = 2;
oo = load 'out_re395_6';
Pdf;
p5 = pdf;
f5 = f;


subplot(2,1,1);
fp = [f1 p1 f2 p2 f3 p3];
semilogx(...
%plot(...
         bins, f1, '1', ...
         bins, p1, '-@01',
         bins, f2, '2', ...
         bins, p2, '-@02',
         bins, f3, '3', ...
         bins, p3, '-@03');
m = max(max([f1 f2 f3 p1 p2 p3]));
axis([0.5*bins(1) 1.5*bins(end) f(1) 1.2*m]);
%axis([0 0.5 0 40]);
title("Vary L, dp/dx");
xlabel("Size");
ylabel("PDF");

subplot(2,1,2);
fp = [f1 p1 f4 p4 f5 p5];
semilogx(...
%plot(...
         bins, f1, '1', ...
         bins, p1, '-@01',
         bins, f4, '2', ...
         bins, p4, '-@02',
         bins, f5, '3', ...
         bins, p5, '-@03');
m = max(max([f1 f4 f5 p1 p4 p5]));
axis([bins(1) bins(end) f(1) 1.2*m]);
%axis([0 0.5 0 10]);
title("Vary mu, dp/dx");
xlabel("Size");
ylabel("PDF");
