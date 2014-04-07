%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guangyuan Sun 05/12
% driver function
% Process homogeneous isotropic turbulence
% Standalone code.
% Mutiple realizations
% Notice make sure you step up correctly in the first block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir = 'homoSL3_040714'; % directory name (or name of data folder)
nrlzn = 5; % # of realization
npart = 3; % # of particles

jrf = 8; % # of particle referece point in dumpTimes.inp (for partDisp.m)

time = 1.6; % evolution time (for partHist.m)
ntime = 1000; % how many interpolation time steps (for partHist.m)
Uinit = 6.55; % initial gas velocity

leftDomain = 0; % left edge of domain (always 0)
rightDomain = 0.508; % right edge of domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partDisp(dir, nrlzn, npart, jrf);
partHist(dir, nrlzn, npart, time, ntime, Uinit);
gasProc(dir, nrlzn, leftDomain, rightDomain, Uinit);

