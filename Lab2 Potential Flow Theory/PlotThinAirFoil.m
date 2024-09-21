function [Press_n, Vel] = PlotThinAirFoil(c,alpha,V_inf,p_inf,rho_inf,N)
%Plot a thin airfoil 
% Use thin airfoil theory to plot a thin airfoil
%
% Author: Maklen A. Estrada
% Date: October 9th, 2022
%%
%Grid conditions
    Xmax = (3*c)/2;
    Xmin = -c/2;
    Ymax = c/2;
    Ymin = -c/2;
%Change in X&Y
    deltaX = Xmax-Xmin;
    deltaY = Ymax-Ymin;
%Bound
    Xbound = [Xmin Xmax];
    Ybound = [Ymin Ymax];
    bound = [Xbound Ybound];
%Step size
    Xn = 100;
    Yn = 100;
%Grid
    [X,Y] = meshgrid(linspace(Xmin,Xmax,Xn),linspace(Ymin,Ymax,Yn));
    
%% Circulation Calculation 
%Seperation Distance
    dX = c./N;
    Xvortex = linspace(dX/2,c-dX,N);
%Strength
    gamma = 2*alpha*V_inf*sqrt( (1-(Xvortex/c))./(Xvortex/c) );
    Circulation = gamma.*dX;
%Radius
rad = @(x1) sqrt((X-x1).^2+(Y).^2);
%% Phi and Psi
%VP for Vortex Flow and SL for Vortex flow
Phi_v = 0;
Psi_v = 0;
for i = 1:N
    theta = atan2(-Y,-X+Xvortex(i));
    Phi_v = Phi_v - (Circulation(i)*theta)/(2*pi);% add thing here
    Psi_v = Psi_v + (Circulation(i)*log(rad(Xvortex(i))))/(2*pi);
end
%VP for uniform flow
Phi_u = (X*cos(alpha) - Y*sin(alpha))*V_inf;
%Vortex and Uniform Flow
Phi = Phi_v + Phi_u;
%SL for uniform flow
Psi_u = (Y*cos(alpha) - X*sin(alpha))*V_inf;
%Vortex and Uniform Flow
Psi = Psi_v + Psi_u;
%% Pressure
%Dynamic Pressure
q = (rho_inf*(V_inf^2))/2;
%Pressure Calculation
Cp = 1 - (gradient(Phi, deltaX/Xn, (deltaY/Yn)./V_inf)).^2;
Press = p_inf + Cp + q;
Press_n = norm(Press);
Cp_n = norm(Cp);
Vel = (1-sqrt(Cp_n))*V_inf;
%% Plots 
%Streamline Plot
figure
subplot(1,3,1)
contourf(X,Y,Psi,90,'k--')
hold on
plot([0 c],[0 0],'k','linewidth',3)
axis(bound)
xlabel('X')
ylabel('Y')
title(['Stream Lines for ' num2str(N) ' vorticies'])
hold off 

%Equipotential Plot
subplot(1,3,2)
contourf(X,Y,Phi,90,'k--')
%contourf(X,Y,Phi,80)
hold on 
plot([0 c],[0 0],'k','linewidth',3)
axis(bound)
xlabel('X')
ylabel('Y')
title(['Equipotential Lines for ' num2str(N) ' vorticies'])
hold off

%Pressure Contour
subplot(1,3,3)
contourf(X,Y,Press,90,'k--')
hold on 
plot([0 c],[0 0],'k','linewidth',3)
axis(bound)
axis(bound)
xlabel('X')
ylabel('Y')
title(['Pressure Contour Lines for ' num2str(N) ' vorticies'])
hold off

set(gcf,'Position', [584,432,1741,583])