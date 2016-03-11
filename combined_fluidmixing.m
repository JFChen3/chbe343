function combined_fluidmixing
close all

%obtain vals
[Re_6265, P_6265, Re_6243, P_6243, Re_6221, P_6221] = experiment();
[Re, Povec1, Povec2, Povec3, ~, ~, ~, ~] = po_derivation();

%plot experimental data to derived Po
figure
loglog(Re_6265, P_6265, 'ro-')
hold on
loglog(Re_6243, P_6243, 'g^-')
loglog(Re_6221, P_6221, 'b.-')
loglog(Re, Povec1, 'r--')
loglog(Re, Povec2, 'g--')
loglog(Re, Povec3, 'b--')
axis([200 3000 0.8 11])
legend('62x65 mm blade','62x43 mm blade','62x21 mm blade', 'Np for 62x65 mm blade', 'Np for 62x43 mm blade', 'Np for 62x21 mm blade', 'Location', 'Northeast');
ylabel('Power Number')
xlabel('Reynolds Number')

end

function [Re_6265, P_6265, Re_6243, P_6243, Re_6221, P_6221] = experiment()

%Given lever arm distance (m)
lever_arm = 0.11;

% Rotational speed in rpm, force in N
% Impeller size 62x65 mm
rpm_6265 = [64.5, 114.5, 134.0, 152.1, 170.5, 195.0, 212.0]; 
force_6265 = [0.45, 2.0, 3.8, 5.3, 6.0, 7.9, 9.9];
[rot_speed_6265, torque_6265, power_6265] = calc_params(rpm_6265, force_6265, lever_arm);
disp('---------------------------------------------------------------------------')
disp('Calculated Parameters, 62x65mm impeller')
disp('---------------------------------------------------------------------------')
disp('Rotation, rpm   Force, N   Rotational Speed, rad/s   Torque, N*m   Power, W')
for i = 1:numel(rpm_6265)
    fprintf('%.1f             %.2f       %.1f                      %.2f          %.1f \n',rpm_6265(i), force_6265(i), rot_speed_6265(i), torque_6265(i), power_6265(i))
end
disp('---------------------------------------------------------------------------')

% Impeller size 62x43 mm
rpm_6243 = [46.3, 95.4, 130.7, 153.9, 174.9, 199.5, 235.4, 250.0];
force_6243 = [1.1, 1.9, 3.1, 4.2, 4.8, 6.3, 8.6, 9.6];
[rot_speed_6243, torque_6243, power_6243] = calc_params(rpm_6243, force_6243, lever_arm);
disp('---------------------------------------------------------------------------')
disp('Calculated Parameters, 62x43mm impeller')
disp('---------------------------------------------------------------------------')
disp('Rotation, rpm   Force, N   Rotational Speed, rad/s   Torque, N*m   Power, W')
for i = 1:numel(rpm_6243)
    fprintf('%.1f             %.2f       %.1f                      %.2f          %.1f \n',rpm_6243(i), force_6243(i), rot_speed_6243(i), torque_6243(i), power_6243(i))
end
disp('---------------------------------------------------------------------------')

% Impeller size 62x21 mm
rpm_6221 = [72.0, 93.1, 128.0, 158.5, 183.0, 208.0, 252.0, 287.0, 309.4, 338.0, 367.1, 390.3];
force_6221 = [0.7, 0.7, 1.7, 2.4, 2.6, 3.4, 4.7, 5.3, 6.3, 7.4, 8.5, 9.4];
[rot_speed_6221, torque_6221, power_6221] = calc_params(rpm_6221, force_6221, lever_arm);
disp('---------------------------------------------------------------------------')
disp('Calculated Parameters, 62x21mm impeller')
disp('---------------------------------------------------------------------------')
disp('Rotation, rpm   Force, N   Rotational Speed, rad/s   Torque, N*m   Power, W')
for i = 1:numel(rpm_6221)
    fprintf('%.1f             %.2f       %.1f                      %.2f          %.1f \n',rpm_6221(i), force_6221(i), rot_speed_6221(i), torque_6221(i), power_6221(i))
end
disp('---------------------------------------------------------------------------')

%Power vs speed plot
figure
hold on
plot(rot_speed_6265, power_6265, 'ro-', rot_speed_6243, power_6243, 'g^-', rot_speed_6221, power_6221, 'b.-')
xlabel('Rotational speed, rad/s')
ylabel('Power, W')
h = legend('62 x 65','62 x 43','62 x 21','Location','Southeast');
v = get(h,'title');
set(v,'string','Impeller Size, mm');

%Power Number vs Reynolds Number
density = 865; %mineral oil, kg/m^3
viscosity = 110e-6*density; %Pa*s
% all_rps = 60*[rpm_6265 rpm_6243 rpm_6221];
% all_power = [power_6265 power_6243 power_6221];
[Re_6265, P_6265] = dimensionless_groups(0.062*2+0.055, rpm_6265./60, density, viscosity, power_6265);
[Re_6243, P_6243] = dimensionless_groups(0.062*2+0.055, rpm_6243./60, density, viscosity, power_6243);
[Re_6221, P_6221] = dimensionless_groups(0.062*2+0.055, rpm_6221./60, density, viscosity, power_6221);

% [reynolds_number, power_number] = dimensionless_groups(lever_arm*2, all_rps, density, viscosity, all_power);
% reynolds_number
% power_number
figure
loglog(Re_6265, P_6265, 'ro-', Re_6243, P_6243, 'g^-', Re_6221, P_6221, 'b.-')
xlabel('Reynolds Number')
ylabel('Power Number')
axis([200 2000 .8 10])
h = legend('62 x 65','62 x 43','62 x 21','Location','NorthEast');
v = get(h,'title');
set(v,'string','Impeller Size, mm');

end

function [rot_speed, torque, power] = calc_params(rpm, force, lever_arm)
% Calculate rotational speed, torque, and power

rot_speed = 2*pi*rpm./60; %rad/s
torque = force.*lever_arm; %N*m
power = torque.*rot_speed; %W

end

function [reynolds_number, power_number] = dimensionless_groups(D_impeller, N_rps, density, viscosity, power)

reynolds_number = ((D_impeller.^2).*N_rps.*density)./viscosity;
power_number = power./((N_rps.^3) .* (D_impeller.^5) .* density);

end

function [Re, Povec1, Povec2, Povec3, Re_lit1, Re_lit2, Po_lit1, Po_lit2] = po_derivation()

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
% figure
% loglog(Re, Povec1, 'r')
% hold on
% loglog(Re, Povec2, 'g')
% loglog(Re, Povec3, 'b')
% xlabel('Reynolds number')
% ylabel('Power Number')
% axis([0 100000 0.7 11])
% h = legend('62 x 65','62 x 43','62 x 21','Location','NorthEast');
% v = get(h,'title');
% set(v,'string','Impeller Size, mm');


%for literature curves
alpha_lit = 6;
h = 0.01;
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
loglog(Re_lit1, Povec_lit1, 'r')
hold on
loglog(Re_lit2, Povec_lit2, 'g')
xlabel('Reynolds number')
ylabel('Power Number')
axis([10 100000 0.7 11])
h = legend('Curve 2, W/D = 1/5', 'Curve 4, W/D = 1/8','Location','SouthEast');
v = get(h,'title');
set(v,'string','Impeller Geometry');

%Print stuff
disp(' ')
disp('---------------------------------------------------------------------------')
disp('Power Number Derived from Drag Force');
disp('---------------------------------------------------------------------------')
disp('Impeller        Power Number');
fprintf('62 x 65 mm: \t%.2f\n', Po1)
fprintf('62 x 43 mm: \t%.2f\n', Po2)
fprintf('62 x 21 mm: \t%.2f\n', Po3)
fprintf('Curve 2:\t\t%.3f\n', Po_lit1)
fprintf('Curve 4:\t\t%.3f\n', Po_lit2)
disp(' ')

end

function [Re] = reynolds(D, N, rho, mu)

Re = (D.^2).*N.*rho./mu;
end

function Po = powernum(alpha, Cd, h, r2, r1)

numer = alpha*Cd*h*((r2^4)-(r1^4))*(pi^3);
denom = (2*r2)^5;
Po = numer./denom;
end