%% ASEN 3111 - Computational Assignment 02 - Main
% In this lab we look at thin and thick airfoils and calculate 
% their aerodynamic forces
% 
%
% Author: Maklen A. Estrada
% Collaborators: Ian Wong, Mathew Davis
% Date: October 9th, 2022
clear all;clc;close all
%% Problem 1

%Constants 
c = 2;
alpha = deg2rad(12);
V_inf = 68;
rho_inf = 1.225;
p_inf = 101.3E3;
N=500;

%Call function to plot thin airfoil
[Pres_500,Vel_500] = PlotThinAirFoil(c,alpha,V_inf,p_inf,rho_inf,N);

%% Problem 2 Conducting a Study
%Varying number of vorticies
N_10 = 10;
[Pres_10, Vel_10] = PlotThinAirFoil(c,alpha,V_inf,p_inf,rho_inf,N_10);
N_100 = 100;
[Pres_100,Vel_100] = PlotThinAirFoil(c,alpha,V_inf,p_inf,rho_inf,N_100);
N_1000 = 1000;
[Pres_1000,Vel_1000] = PlotThinAirFoil(c,alpha,V_inf,p_inf,rho_inf,N_1000);
N_10000 = 10000;
[Pres_10000,Vel_10000]= PlotThinAirFoil(c,alpha,V_inf,p_inf,rho_inf,N_10000);

%Calculating error for pressure
er1 = (Pres_10 - Pres_10000)/Pres_10000*100;
er2 = (Pres_100 - Pres_10000)/Pres_10000*100;
er3 = (Pres_500 - Pres_10000)/Pres_10000*100;
er4 = (Pres_1000 - Pres_10000)/Pres_10000*100;
fprintf('Pressure percent error for 10 vorticies, %f',er1)
fprintf('\n')
fprintf('Pressure percent error for 100 vorticies, %f',er2)
fprintf('\n')
fprintf('Pressure percent error for 500 vorticies, %f',er3)
fprintf('\n')
fprintf('Pressure percent error for 1000 vorticies, %f',er4)
fprintf('\n')

%Calculating error for velocity
er1 = (Vel_10000 - Vel_10)/Vel_10000*100;
er2 = (Vel_10000 - Vel_100)/Vel_10000*100;
er3 = (Vel_10000 - Vel_500)/Vel_10000*100;
er4 = (Vel_10000 - Vel_1000)/Vel_10000*100;
fprintf('Velocity percent error for 10 vorticies, %f',er1)
fprintf('\n')
fprintf('Velocity percent error for 100 vorticies, %f',er2)
fprintf('\n')
fprintf('Velocity percent error for 500 vorticies, %f',er3)
fprintf('\n')
fprintf('Velocity percent error for 1000 vorticies, %f',er4)
fprintf('\n')

fprintf('Figures 2-5 show the effect of 10, 100, 1000, and 10,000 discrete vorticies')
fprintf('\n')

%Making plots to compare number of vorticies to pressure and velocity
Pressure = [Pres_10,Pres_100,Pres_500,Pres_1000];
Vel = [Vel_10,Vel_100,Vel_500,Vel_1000];
Number = [10,100,500,1000];

figure
plot(Number,Pressure)
xlabel('Number of Vorticies')
ylabel('Pressure')
title('Pressure compared to number of vorticies')
yline(Pres_10000,'--','Actual Pressure')
grid on

figure
plot(Number,Vel)
xlabel('Number of Vorticies')
ylabel('Velocity')
title('Velocity compared to number of vorticies')
yline(Vel_10000,'--','Actual Velocity')
grid on

fprintf('Figure 6 shows the pressure compared to the number of vorticies')
fprintf('\n')
fprintf('Figure 7 shows the velocity compared to the number of vorticies')
fprintf('\n')
%% Problem 3.a
%Chord Length variation
chord = [1,1.75,2.5,3.25];
figure
for i = 1:length(chord)
subplot(1,4,i)
StreamLines(chord(i),alpha,V_inf,p_inf,rho_inf,N)
set(gcf,'Position', [205,610,2302,688])
end
figure
for i = 1:length(chord)
subplot(1,4,i)
EquipotentialLines(chord(i),alpha,V_inf,p_inf,rho_inf,N)
set(gcf,'Position', [205,610,2302,688])
end

fprintf('Figures 8-9 show the StreamLines and Equipotential Lines for chords length of 1, 1.75, 2.5, and 3.25 meters respectivly')
fprintf('\n')

%% Problem 3.b
%AoA variation
AoA = deg2rad([9,12,15,18]);
figure
for i = 1:length(AoA)
subplot(1,4,i)
StreamLines(c,AoA(i),V_inf,p_inf,rho_inf,N)
set(gcf,'Position', [205,610,2302,688])
end
figure
for i = 1:length(AoA)
subplot(1,4,i)
EquipotentialLines(c,AoA(i),V_inf,p_inf,rho_inf,N)
set(gcf,'Position', [205,610,2302,688])
end

fprintf('Figures 10-11 show the StreamLines and Equipotential Lines for AoA of 9, 12, 15, and 18 degrees respectivly')
fprintf('\n')

%% Problem 3.c
%Flow Speed variation
V = deg2rad([60,70,80,90]);
figure
for i = 1:length(V)
subplot(1,4,i)
StreamLines(c,alpha,V(i),p_inf,rho_inf,N)
set(gcf,'Position', [205,610,2302,688])
end
figure
for i = 1:length(V)
subplot(1,4,i)
EquipotentialLines(c,alpha,V(i),p_inf,rho_inf,N)
set(gcf,'Position', [205,610,2302,688])
end

fprintf('Figures 12-13 show the StreamLines and Equipotential Lines for Free-Stream Flow Speeds of 60, 70, 80, and 90 m/s respectivly')
fprintf('\n')