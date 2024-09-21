%% ASEN 3111 - Computational Assignment 04 - Main
%In this computational assignment we looked at how AoA and Mach affects the
%sectional wave-drag coefficient and coefficient of lift
%
% Author: Maklen Estrada
% Collaborators: Ian Wong
% Date: December 4th, 2022
clc
clear all
close all

gamma = 1.4;
theta = [0:4:48];
%Only looking at some mach values
Mach = [1,1.5,2,2.5,3];
theta_l = length(theta);
Mach_l = length(Mach);
type_s = 'Strong';
type_w = 'Weak';
for i=1:theta_l
    for j=1:Mach_l
        beta_strong(i,j) = ObliqueShockBeta(Mach(j), theta(i),gamma,type_s);
        beta_weak(i,j) = ObliqueShockBeta(Mach(j), theta(i),gamma,type_w);
    end
end

%Getting the real component of the beta's 
beta_strong = real(beta_strong);
Pos_bs = NaN(size(beta_strong));
idxs = beta_strong > 0;
Pos_bs(idxs) = beta_strong(idxs);

%Makes a new matrix with only positive values
beta_weak = real(beta_weak);
Pos_bw = NaN(size(beta_strong));
idxw = beta_weak > 0;
Pos_bw(idxw) = beta_weak(idxw);

%Theta-Beta-M plot
figure
hold on
for i=2:Mach_l
plot(theta,Pos_bs(:,i),'LineWidth',2);
plot(theta, Pos_bw(:,i),'LineWidth',2);
ylim([0 100]);
end
grid on;
grid minor;
xlabel('\alpha');set(gca,'FontSize',15)
ylabel('\beta');set(gca,'FontSize',15)
title('\theta - \beta - M Diagram')
legend('Mach 1','Mach 1.5','Mach 2','Mach 2.5','Mach 2.5','Mach 3','Location','southeast')
legend('boxoff')
%% Problem 2
[c_l, c_dw] = DiamondAirfoil(2,deg2rad(5),deg2rad(10),deg2rad(5));
fprintf('Probelm 2')
fprintf('\n')
fprintf('Coefficient of Lift: %f',c_l)
fprintf('\n')
fprintf('Sectional Wave-Drag Coefficient: %f',c_dw)
fprintf('\n')
%% Problem 3

M3 = [2:5];
alpha3 = [.2:.2:10];
M3_l = length(M3);
a3_l = length(alpha3);
for i =1:a3_l
    for j=1:M3_l
        [c_l3(i,j), c_dw3(i,j)] = DiamondAirfoil(M3(j),alpha3(i), deg2rad(10), deg2rad(5));
    end
end


%Coefficient of lift plot
figure
hold on
for i = 1:M3_l
plot(alpha3, c_l3(:,i),'LineWidth',2)
xlabel('\alpha');set(gca,'FontSize',15)
ylabel('C_L');set(gca,'FontSize',15)
legend('Mach 2', 'Mach 3','Mach 4', 'Mach 5')
grid on;
end

%Sectional Wave-Drag Coefficient
figure
hold on
for i=1:M3_l
plot(alpha3, c_dw3(:,i),'LineWidth',2)
xlabel('\alpha');set(gca,'FontSize',15);
ylabel('C_{d,w}');set(gca,'FontSize',15);
legend('Mach 2', 'Mach 3','Mach 4', 'Mach 5')
grid on;
end
warning('off','all')