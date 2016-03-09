function fluidmixing
close all

h = 0.065; %m
w = 0.062; %m
r1 = 0.055/2; %m
r2 = r1 + w; %m
alpha = 2; %# of blades
Cd = 2; %from drag coeffs doc
rho = 864; %kg/m^3
mu = 0.0952; %kg/m*s

N = 0:250; %rev/sec
Re = reynolds(r2, N, rho, mu);
Po = powernum(alpha, Cd, h, r2, r1, N, rho, Re, mu);

%the plot
loglog(Re, Po)

end

function [Re] = reynolds(r2, N, rho, mu)

Re = ((2*r2).^2).*N.*rho./mu;
end

function Po = powernum(alpha, Cd, h, r2, r1, N, rho, Re, mu)

numer = alpha*Cd*h*((r2^3)-(r1^3))*(r2-r1);
denom = 48*(pi^3)*((2*r2)^5);

Po = numer./denom;
end
