function fluidmixing
close all

%constants
alpha = 2; %# of blades
Cd = 2; %from drag coeffs doc
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
axis([0 100000 0.7 11])
h = legend('62 x 65','62 x 43','62 x 21','Location','NorthEast');
v = get(h,'title');
set(v,'string','Impeller Size, mm');

%Print stuff
disp('             Power Number');
fprintf('62 x 65 mm: \t%.2f\n', Po1)
fprintf('62 x 43 mm: \t%.2f\n', Po2)
fprintf('62 x 21 mm: \t%.2f\n', Po3)

%for literature curves
alpha_lit = 6;
h = 1;
r1 = 0;
D_lit1 = 5*h;
r2_lit1 = D_lit1/2;
D_lit2 = 8*h;
r2_lit2 = D_lit2/2;

%curve 2, W/D = 1/5
Re_lit1 = reynolds(D_lit1, N, rho, mu);
Po_lit1 = powernum(alpha_lit, Cd, h, r2_lit1, r1);
Povec_lit1 = Po_lit1*ones(size(Re_lit1));

%curve 4, W/D = 1/8
Re_lit2 = reynolds(D_lit2, N, rho, mu);
Po_lit2 = powernum(alpha_lit, Cd, h, r2_lit2, r1);
Povec_lit2 = Po_lit2*ones(size(Re_lit2));

%the plotting
figure
loglog(Re, Povec_lit1, 'r')
hold on
loglog(Re, Povec_lit2, 'g')
xlabel('Reynolds number')
ylabel('Power Number')
axis([0 100000 0.7 11])
h = legend('Curve 2, W/D = 1/5', 'Curve 4, W/D = 1/8','Location','SouthEast');
v = get(h,'title');
set(v,'string','Impeller Geometry');

fprintf('Curve 2: \t%.2f\n', Po_lit1)
fprintf('Curve 3: \t%.2f\n', Po_lit2)

% For the 8 bladed impeller
alpha_8 = 8;
r1_8 = .051;
r2_8 = .071;
h_8 = .012;

omega_8 = 0:45;
N_8 = omega_8 ./ (2 * pi);
D_8 = 2 * r2_8;

Po_8 = powernum(alpha_8, Cd, h_8, r2_8, r1_8);
Power_8 = Po_8 .* rho .* N_8.^3 .* D_8.^5;

Re_8 = reynolds(D_8, N_8, rho, mu);

% Plotting
figure()
hold on
plot(omega_8,Power_8)
xlabel('Rotational Speed, rad/s')
ylabel('Power, W')

figure()
Re_Vec = 1e3:1e5;
Po_Vec = ones(1,numel(Re_Vec)).*(Po_8);
plot(Re_Vec, Po_Vec)
xlabel('Reynolds Number')
ylabel('Power Number')

end

function [Re] = reynolds(D, N, rho, mu)

Re = (D.^2).*N.*rho./mu;
end

function Po = powernum(alpha, Cd, h, r2, r1)

numer = alpha*Cd*h*((r2^4)-(r1^4))*(pi^3);
denom = (2*r2)^5;
Po = numer./denom;
end
