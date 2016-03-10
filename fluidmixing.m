function fluidmixing
close all

%constants
alpha = 2; %# of blades
Cd = 1.28; %from drag coeffs doc
rho = 864; %kg/m^3
mu = 0.0952; %kg/m*s

%blade dimensions
h1 = 0.065; %m
h2 = 0.043; %m
h3 = 0.021; %m
w = 0.062; %m
r1 = 0.055/2; %m
r2 = r1 + w; %m
D = 2*r2;

%for literature
% r1 = 0;
% D = 8*h;
% r2 = D/2;
% w = r2 - r1;

N = linspace(0.1, 10^5); %rev/sec
Re = reynolds(D, N, rho, mu);

Po1 = powernum(alpha, Cd, h1, r2, r1);
Povec1 = Po1*ones(size(Re));

Po2 = powernum(alpha, Cd, h2, r2, r1);
Povec2 = Po2*ones(size(Re));

Po3 = powernum(alpha, Cd, h3, r2, r1);
Povec3 = Po3*ones(size(Re));

%the plotting
figure
loglog(Re, Povec1, 'r')
hold on
loglog(Re, Povec2, 'g')
loglog(Re, Povec3, 'b')
xlabel('Reynolds number')
ylabel('Power Number')
title('Power Number Derived from Drag Force')
axis([0 100000 0.4 5])
h = legend('62 x 65','62 x 43','62 x 21','Location','NorthEast');
v = get(h,'title');
set(v,'string','Impeller Size, mm');

end

function [Re] = reynolds(D, N, rho, mu)

Re = (D.^2).*N.*rho./mu;
end

function Po = powernum(alpha, Cd, h, r2, r1)

numer = alpha*Cd*h*((r2^4)-(r1^4))*(pi^3);
denom = (2*r2)^5;
Po = numer./denom;
end
