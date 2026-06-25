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
h = 0.001; %Integration Step
r01 = 120; %mm
r02 = 80; %mm
r03 = 100; %mm
r04 = 150; %mm
W1 = 10; %Hz
W2 = 15; %Hz
W3 = 8; %Hz
W4 = 5; %Hz

%Load Optistruct Matrices
load 'Frame_M.dat'
M1 = Frame_M;
load 'Frame_K.dat'
K1 = Frame_K;
load 'Frame_C.dat'
C1 = Frame_C;

%Construct Global Mass, Stiffness and Damping Matrices
dof = 66 + 4; %Total Dof
M = zeros(dof); %Global Mass Matrix
K = zeros(dof); %Global Stiffness Matrix
C = zeros(dof); %Global Damping Matrix
M(1:66,1:66) = M(1:66,1:66) + M1; %Mass Matrix from the Basic Chassis Structure
%Mass from the Suspension-Wheel Subsystems
M(67,67) = M(67,67) + mw;
M(68,68) = M(68,68) + mw;
M(69,69) = M(69,69) + mw;
M(70,70) = M(70,70) + mw;
K(1:66,1:66) = K(1:66,1:66) + K1; %Stiffness Matrix from the Basic Chassis Structure
%Stiffness from the Suspension-Wheel Subsystems
K([63,67],[63,67]) = K([63,67],[63,67]) + [ks -ks; -ks ks+kw];
K([66,68],[66,68]) = K([66,68],[66,68]) + [ks -ks; -ks ks+kw];
K([33,69],[33,69]) = K([33,69],[33,69]) + [ks -ks; -ks ks+kw];
K([36,70],[36,70]) = K([36,70],[36,70]) + [ks -ks; -ks ks+kw];
C(1:66,1:66) = C(1:66,1:66) + C1; %Damping Matrix from the Basic Chassis Structure
%Damping from the Suspension-Wheel Subsystems
C([63,67],[63,67]) = C([63,67],[63,67]) + [cs -cs; -cs cs+cw];
C([66,68],[66,68]) = C([66,68],[66,68]) + [cs -cs; -cs cs+cw];
C([33,69],[33,69]) = C([33,69],[33,69]) + [cs -cs; -cs cs+cw];
C([36,70],[36,70]) = C([36,70],[36,70]) + [cs -cs; -cs cs+cw];

%Solve Transient Problem with Numerical Integration Method 'Central Difference'
u = zeros(2003,70); %Total Displacement Matrix
v = zeros(2001,70); %Total Velocity Matrix
a = zeros(2001,70); %Total Acceleration Matrix
%Ground Excitation and the Respective Derivative on each Suspension-Wheel Subsystem at Time Instance 0
R1_0 = r01*cos(2*pi*W1*0);
R1d_0 = -2*pi*W1*r01*sin(2*pi*W1*0);
R2_0 = r02*sin(2*pi*W2*0);
R2d_0 = 2*pi*W2*r02*cos(2*pi*W2*0);
R3_0 = r03*sin(2*pi*W3*0);
R3d_0 = 2*pi*W3*r03*cos(2*pi*W3*0);
R4_0 = r04*cos(2*pi*W4*0);
R4d_0 = -2*pi*W4*r04*sin(2*pi*W4*0);
%Excitation Force on each Suspension-Wheel Subsystem at Time Instance 0
F_0 = zeros(70,1);
F_0(67) = cw*R1d_0 + kw*R1_0;
F_0(68) = cw*R2d_0 + kw*R2_0;
F_0(69) = cw*R3d_0 + kw*R3_0;
F_0(70) = cw*R4d_0 + kw*R4_0;
a(1,:) = (inv(M)*(F_0 - C*v(1,:)' - K*u(2,:)'))'; %Acceleration at Time Instance 0
u(1,:) = 1/2*(h^2*a(1,:) - 2*h*v(1,:) + 2*u(2,:)); %Displacement at Time Instance -1
t = h;
A = inv(M/h^2 + C/2/h + K/3); %Compute Stable Matrix A for Displacement Calculation at Time Instance t
while t <= 2+h
    %Ground Excitation and the Respective Derivative on each Suspension-Wheel Subsystem at Time Instance t-2h
    R1_1 = r01*cos(2*pi*W1*(t-2*h));
    R1d_1 = -2*pi*W1*r01*sin(2*pi*W1*(t-2*h)); 
    R2_1 = r02*sin(2*pi*W2*(t-2*h));
    R2d_1 = 2*pi*W2*r02*cos(2*pi*W2*(t-2*h)); 
    R3_1 = r03*sin(2*pi*W3*(t-2*h));
    R3d_1 = 2*pi*W3*r03*cos(2*pi*W3*(t-2*h));
    R4_1 = r04*cos(2*pi*W4*(t-2*h));
    R4d_1 = -2*pi*W4*r04*sin(2*pi*W4*(t-2*h));
    %Excitation Force on each Suspension-Wheel Subsystem at Time Instance t-2h
    F_1 = zeros(70,1); 
    F_1(67) = cw*R1d_1 + kw*R1_1;
    F_1(68) = cw*R2d_1 + kw*R2_1;
    F_1(69) = cw*R3d_1 + kw*R3_1;
    F_1(70) = cw*R4d_1 + kw*R4_1;
    %Ground Excitation and the Respective Derivative on each Suspension-Wheel Subsystem at Time Instance t-h
    R1_2 = r01*cos(2*pi*W1*(t-h));
    R1d_2 = -2*pi*W1*r01*sin(2*pi*W1*(t-h));
    R2_2 = r02*sin(2*pi*W2*(t-h));
    R2d_2 = 2*pi*W2*r02*cos(2*pi*W2*(t-h));
    R3_2 = r03*sin(2*pi*W3*(t-h));
    R3d_2 = 2*pi*W3*r03*cos(2*pi*W3*(t-h));
    R4_2 = r04*cos(2*pi*W4*(t-h));
    R4d_2 = -2*pi*W4*r04*sin(2*pi*W4*(t-h));
    %Excitation Force on each Suspension-Wheel Subsystem at Time Instance t-h
    F_2 = zeros(70,1);
    F_2(67) = cw*R1d_2 + kw*R1_2;
    F_2(68) = cw*R2d_2 + kw*R2_2;
    F_2(69) = cw*R3d_2 + kw*R3_2;
    F_2(70) = cw*R4d_2 + kw*R4_2;
    %Ground Excitation and the Respective Derivative on each Suspension-Wheel Subsystem at Time Instance t
    R1_3 = r01*cos(2*pi*W1*t);
    R1d_3 = -2*pi*W1*r01*sin(2*pi*W1*t);
    R2_3 = r02*sin(2*pi*W2*t);
    R2d_3 = 2*pi*W2*r02*cos(2*pi*W2*t);
    R3_3 = r03*sin(2*pi*W3*t);
    R3d_3 = 2*pi*W3*r03*cos(2*pi*W3*t);
    R4_3 = r04*cos(2*pi*W4*t);
    R4d_3 = -2*pi*W4*r04*sin(2*pi*W4*t);
    %Excitation Force on each Suspension-Wheel Subsystem at Time Instance t
    F_3 = zeros(70,1);
    F_3(67) = cw*R1d_3 + kw*R1_3;
    F_3(68) = cw*R2d_3 + kw*R2_3;
    F_3(69) = cw*R3d_3 + kw*R3_3;
    F_3(70) = cw*R4d_3 + kw*R4_3;

    u(round(t/h)+2,:) = (A*(1/3*(F_3+F_2+F_1) + (2*M/h^2-K/3)*u(round(t/h)+1,:)' + (-M/h^2+C/2/h-K/3)*u(round(t/h),:)'))'; %Displacement

    if t/h >= 2
        v(round(t/h),:) = 1/2/h*(u(round(t/h)+2,:) - u(round(t/h),:)); %Velocity
        a(round(t/h),:) = 1/h^2*(u(round(t/h)+2,:) - 2*u(round(t/h)+1,:) + u(round(t/h),:)); %Acceleration
    end

    t = t + h;
end

%Comparison Diagrams for Point Β1
load 'export_data_B_1.col'
B_1_plot = export_data_B_1;
t= 0:0.001:2;

f1 = figure;
hold on
plot(t, u(2:2002,63), "r");
plot(B_1_plot(1:2001,1), B_1_plot(1:2001,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of B1 on z axis");
grid on
hold off
f2 = figure;
hold on
plot(t, v(:,63), "r");
plot(B_1_plot(1:2001,1), B_1_plot(2002:4002,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of B1 on z axis");
grid on
hold off
f3 = figure;
hold on
plot(t, a(:,63), "r");
plot(B_1_plot(1:2001,1), B_1_plot(4003:6003,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of B1 on z axis");
grid on
hold off

%Comparison Diagrams for Point Β2
load 'export_data_B_2.col'
B_2_plot = export_data_B_2;

f4 = figure;
hold on
plot(t, u(2:2002,66), "r");
plot(B_2_plot(1:2001,1), B_2_plot(1:2001,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of B2 on z axis");
grid on
hold off
f5 = figure;
hold on
plot(t, v(:,66), "r");
plot(B_2_plot(1:2001,1), B_2_plot(2002:4002,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of B2 on z axis");
grid on
hold off
f6 = figure;
hold on
plot(t, a(:,66), "r");
plot(B_2_plot(1:2001,1), B_2_plot(4003:6003,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of B2 on z axis");
grid on
hold off

%Comparison Diagrams for Point Β3
load 'export_data_B_3.col'
B_3_plot = export_data_B_3;

f7 = figure;
hold on
plot(t, u(2:2002,33), "r");
plot(B_3_plot(1:2001,1), B_3_plot(1:2001,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of B3 on z axis");
grid on
hold off
f8 = figure;
hold on
plot(t, v(:,33), "r");
plot(B_3_plot(1:2001,1), B_3_plot(2002:4002,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of B3 on z axis");
grid on
hold off
f9 = figure;
hold on
plot(t, a(:,33), "r");
plot(B_3_plot(1:2001,1), B_3_plot(4003:6003,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of B3 on z axis");
grid on
hold off

%Comparison Diagrams for Point Β4
load 'export_data_B_4.col'
B_4_plot = export_data_B_4;

f10 = figure;
hold on
plot(t, u(2:2002,36), "r");
plot(B_4_plot(1:2001,1), B_4_plot(1:2001,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of B4 on z axis");
grid on
hold off
f11 = figure;
hold on
plot(t, v(:,36), "r");
plot(B_4_plot(1:2001,1), B_4_plot(2002:4002,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of B4 on z axis");
grid on
hold off
f12 = figure;
hold on
plot(t, a(:,36), "r");
plot(B_4_plot(1:2001,1), B_4_plot(4003:6003,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of B4 on z axis");
grid on
hold off

%Comparison Diagrams for Point Ι1
load 'export_data_I_1.col'
I_1_plot = export_data_I_1;

f13 = figure;
hold on
plot(t, u(2:2002,55), "r");
plot(I_1_plot(1:2001,1), I_1_plot(1:2001,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I1 on x axis");
grid on
hold off
f14 = figure;
hold on
plot(t, v(:,55), "r");
plot(I_1_plot(1:2001,1), I_1_plot(2002:4002,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I1 on x axis");
grid on
hold off
f15 = figure;
hold on
plot(t, a(:,55), "r");
plot(I_1_plot(1:2001,1), I_1_plot(4003:6003,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I1 on x axis");
grid on
hold off
f16 = figure;
hold on
plot(t, u(2:2002,56), "r");
plot(I_1_plot(1:2001,1), I_1_plot(6004:8004,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I1 on y axis");
grid on
hold off
f17 = figure;
hold on
plot(t, v(:,56), "r");
plot(I_1_plot(1:2001,1), I_1_plot(8005:10005,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I1 on y axis");
grid on
hold off
f18 = figure;
hold on
plot(t, a(:,56), "r");
plot(I_1_plot(1:2001,1), I_1_plot(10006:12006,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I1 on y axis");
grid on
hold off
f19 = figure;
hold on
plot(t, u(2:2002,57), "r");
plot(I_1_plot(1:2001,1), I_1_plot(12007:14007,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I1 on z axis");
grid on
hold off
f20 = figure;
hold on
plot(t, v(:,57), "r");
plot(I_1_plot(1:2001,1), I_1_plot(14008:16008,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I1 on z axis");
grid on
hold off
f21 = figure;
hold on
plot(t, a(:,57), "r");
plot(I_1_plot(1:2001,1), I_1_plot(16009:18009,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I1 on z axis");
grid on
hold off

%Comparison Diagrams for Point Ι2
load 'export_data_I_2.col'
I_2_plot = export_data_I_2;

f22 = figure;
hold on
plot(t, u(2:2002,43), "r");
plot(I_2_plot(1:2001,1), I_2_plot(1:2001,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I2 on x axis");
grid on
hold off
f23 = figure;
hold on
plot(t, v(:,43), "r");
plot(I_2_plot(1:2001,1), I_2_plot(2002:4002,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I2 on x axis");
grid on
hold off
f24 = figure;
hold on
plot(t, a(:,43), "r");
plot(I_2_plot(1:2001,1), I_2_plot(4003:6003,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I2 on x axis");
grid on
hold off
f25 = figure;
hold on
plot(t, u(2:2002,44), "r");
plot(I_2_plot(1:2001,1), I_2_plot(6004:8004,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I2 on y axis");
grid on
hold off
f26 = figure;
hold on
plot(t, v(:,44), "r");
plot(I_2_plot(1:2001,1), I_2_plot(8005:10005,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I2 on y axis");
grid on
hold off
f27 = figure;
hold on
plot(t, a(:,44), "r");
plot(I_2_plot(1:2001,1), I_2_plot(10006:12006,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I2 on y axis");
grid on
hold off
f28 = figure;
hold on
plot(t, u(2:2002,45), "r");
plot(I_2_plot(1:2001,1), I_2_plot(12007:14007,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I2 on z axis");
grid on
hold off
f29 = figure;
hold on
plot(t, v(:,45), "r");
plot(I_2_plot(1:2001,1), I_2_plot(14008:16008,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I2 on z axis");
grid on
hold off
f30 = figure;
hold on
plot(t, a(:,45), "r");
plot(I_2_plot(1:2001,1), I_2_plot(16009:18009,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I2 on z axis");
grid on
hold off

%Comparison Diagrams for Point Ι3
load 'export_data_I_3.col'
I_3_plot = export_data_I_3;

f31 = figure;
hold on
plot(t, u(2:2002,37), "r");
plot(I_3_plot(1:2001,1), I_3_plot(1:2001,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I3 on x axis");
grid on
hold off
f32 = figure;
hold on
plot(t, v(:,37), "r");
plot(I_3_plot(1:2001,1), I_3_plot(2002:4002,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I3 on x axis");
grid on
hold off
f33 = figure;
hold on
plot(t, a(:,37), "r");
plot(I_3_plot(1:2001,1), I_3_plot(4003:6003,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I3 on x axis");
grid on
hold off
f34 = figure;
hold on
plot(t, u(2:2002,38), "r");
plot(I_3_plot(1:2001,1), I_3_plot(6004:8004,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I3 on y axis");
grid on
hold off
f35 = figure;
hold on
plot(t, v(:,38), "r");
plot(I_3_plot(1:2001,1), I_3_plot(8005:10005,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I3 on y axis");
grid on
hold off
f36 = figure;
hold on
plot(t, a(:,38), "r");
plot(I_3_plot(1:2001,1), I_3_plot(10006:12006,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I3 on y axis");
grid on
hold off
f37 = figure;
hold on
plot(t, u(2:2002,39), "r");
plot(I_3_plot(1:2001,1), I_3_plot(12007:14007,2), "b");
xlabel("Time (s)");
ylabel("Displacement (mm)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for displacement of I3 on z axis");
grid on
hold off
f38 = figure;
hold on
plot(t, v(:,39), "r");
plot(I_3_plot(1:2001,1), I_3_plot(14008:16008,2), "b");
xlabel("Time (s)");
ylabel("Velocity (mm/s)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for velocity of I3 on z axis");
grid on
hold off
f39 = figure;
hold on
plot(t, a(:,39), "r");
plot(I_3_plot(1:2001,1), I_3_plot(16009:18009,2), "b");
xlabel("Time (s)");
ylabel("Acceleration (mm/s^2)");
legend(["Central Difference", "Optistruct"]);
title("Comparison for acceleration of I3 on z axis");
grid on
hold off
