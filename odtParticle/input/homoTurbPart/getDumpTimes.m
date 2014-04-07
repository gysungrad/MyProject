clc; clear;

time = 1.6;   % time domain
N = 41;         % number of dump times


timestep = linspace(0,time,N);

outfilename = 'dumpTimes.inp';

%open the output file for write
[fiout,message]=fopen(outfilename,'w');
if fiout<0
    disp(message);
    exit;
end

%output the title
fprintf(fiout, '%d       number of dump times - times follow\n',N);
for j = 1:N
    fprintf(fiout,'%6.10e\n',timestep(j));
end

fclose(fiout);

