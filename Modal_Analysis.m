clear
clc
format long

%Add shared function folder to path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'data')));

%Define Constants
mw = 15.2*10^(-3); %tn
ks = 23.5; %N/mm
kw = 40.5; %N/mm
cs = 0.35; %Ns/mm
cw = 0.15; %Ns/mm

%Load Optistruct Matrices
load 'Frame_M.dat'
M1 = Frame_M;
load 'Frame_K.dat'
K1 = Frame_K;
load 'Frame_C.dat'
C1 = Frame_C;

%Construct Total Mass, Stiffness and Damping Matrices
dof = 66 + 4; %Total Dof
M = zeros(dof); %Total Mass Matrix
K = zeros(dof); %Total Stiffness Matrix
C = zeros(dof); %Total Damping Matrix
M(1:66,1:66) = M(1:66,1:66) + M1; %Mass Matrix from the basic structure 
M(67,67) = M(67,67) + mw; %Mass from the First Subsystem
M(68,68) = M(68,68) + mw; %Mass from the Second Subsystem
M(69,69) = M(69,69) + mw; %Mass from the Third Subsystem
M(70,70) = M(70,70) + mw; %Mass from the Fourth Subsystem
K(1:66,1:66) = K(1:66,1:66) + K1; %Stiffness Matrix from the Basic Structure
K([63,67],[63,67]) = K([63,67],[63,67]) + [ks -ks; -ks ks+kw]; %Stiffness from the First Subsystem
K([66,68],[66,68]) = K([66,68],[66,68]) + [ks -ks; -ks ks+kw]; %Stiffness from the Second Subsystem
K([33,69],[33,69]) = K([33,69],[33,69]) + [ks -ks; -ks ks+kw]; %Stiffness from the Third Subsystem
K([36,70],[36,70]) = K([36,70],[36,70]) + [ks -ks; -ks ks+kw]; %Stiffness from the Fourth Subsystem
C(1:66,1:66) = C(1:66,1:66) + C1; %Damping Matrix from the Basic Structure
C([63,67],[63,67]) = C([63,67],[63,67]) + [cs -cs; -cs cs+cw]; %Damping from the First Subsystem
C([66,68],[66,68]) = C([66,68],[66,68]) + [cs -cs; -cs cs+cw]; %Damping from the Second Subsystem
C([33,69],[33,69]) = C([33,69],[33,69]) + [cs -cs; -cs cs+cw]; %Damping from the Third Subsystem
C([36,70],[36,70]) = C([36,70],[36,70]) + [cs -cs; -cs cs+cw]; %Damping from the Fourth Subsystem

%Eigenproblem
[eigv, wa] = eig(K, M);
[eigw,indx] = sort(sqrt(diag(wa)));
eigv = eigv(:,indx);
disp("First 15 natural frequencies computed via Matlab are (in Hz):");
disp(eigw(1:15)/2/pi);

%Comparison with Optistruct Results (the third natural frequency is omitted because it is a complex number due to the asymmetric stiffness matrix loaded from Optistruct)
opt = [1.00862*10^(-3); 1.064405*10^(-3); 2.736431; 3.539792; 3.703959; 10.2976; 10.51411; 10.74477; 10.78271; 25.18493; 42.48746; 43.39948; 49.93391; 60.53719];
mat = [eigw(1:2)/2/pi; eigw(4:15)/2/pi];
dif = mat - opt;
disp("Difference between natural frequencies computed via Matlab and Optistruct (in Hz):");
disp(dif);
