function StreamLines(c,alpha,V_inf,p_inf,rho_inf,N)
%Plot a thin airfoil 
% Use thin airfoil theory to plot a Stream Lines
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

%% Plots 
%Streamline Plot
contourf(X,Y,Psi,90,'k--')
hold on
plot([0 c],[0 0],'k','linewidth',c)
axis(bound)
xlabel('X')
ylabel('Y')
title(['Stream Lines'])
hold off 