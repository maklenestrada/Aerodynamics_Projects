function [x,y] = NACA(m,p,t,c,N)
%This function calculates the lift slope and zero-lift agnle of
%attack
%
%Author: Maklen A. Estrada
%Collaborators: Ian Wong
%Date: November 6th, 2022

%Thickness distribution of the airfoil
x_t = linspace(0,c,N);
y_t = (t/.2)*c .* ( .2969*sqrt(x_t/c) - .1260*(x_t/c) - .3516*((x_t/c).^2) + .2843*((x_t/2).^3) - .1036*((x_t/c).^4) );

%Formula for the mean camber line
if (p == 0) && (m==0)
    y_c = zeros(1,length(x_t));
else
for i = 1:length(x_t)
if x_t(i) <= p*c
    y_c(i) = m.*(x_t(i)/(p^2))*(2*p - (x_t(i)/c));
elseif x_t(i) >= (p*c)
    y_c(i) = m*( (c - x_t(i))/((1-p)^2) )*(1 + (x_t(i)/c) - 2*p);
end
end
end

%Camber line derrivative
if (p == 0) && (m==0)
    dy_c = zeros(1,length(x_t));
else
for i = 1:length(x_t)
if x_t(i) <= p*c
    dy_c(i) = 2*m .* ( (p*c - x_t(i) ./ (p^2*c)));
elseif x_t(i) >= (p*c)
    dy_c(i) = 2*m .* ( (p*c - x_t(i)) ./ ((1-p)^2*c) );
end
end
end

xi = atan(dy_c);

%Upper Airfoil surface
x_upp = x_t - y_t.*sin(xi);
y_upp = y_c + y_t.*cos(xi);

%Lower Airfoil surface
x_low = x_t + y_t.*sin(xi);
y_low = y_c - y_t.*cos(xi);

%Total Airfoil surface
x = [flip(x_low),x_upp(2:end)];
y = [flip(y_low),y_upp(2:end)];

end