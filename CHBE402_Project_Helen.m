function CHBE402_Project_Helen
close all

%constants
g = 9.81; %m/s^2
R = 8.314; %J/molK
Tm = 166 + 273.15; %MP of IL
dHfus = -19900; %J/mol
dH = -58500; %J/mol
dS = -143; %J/mol
MW_CO2 = 44; %g/mol
D = 0.002; %m

m_in = 596000; %kg/hr
n_in = m_in*1000/MW_CO2; %mol/hr

% Initialize everything
xn = 0.75;
T = 273.15 + 125; %K
error = 100;
tol = 0.001;
dT = 0.01; %K, temp step
while error > tol
    T = T + dT;
    [K, ~, Pcrit] = solve_Pcrit(T, R, Tm, dHfus, dH, dS);
    sol = co2_solubility((T-273.15), Pcrit);
    xntemp = sol/K;
    error = abs(xn-xntemp)/xn;
end

fprintf('Temperature: %.2f K\n', T)

rho_IL = IL_density(T);

%Open file for CO2 data: viscosity and density
%Based on T = 415.27
CO2_fileID = fopen('CO2_data.txt');
CO2_data = fscanf(CO2_fileID, '%f', [13, 91])';

P_CO2 = CO2_data(:, 2)'; %bar
mus_CO2 = CO2_data(:, 12)'; %Pa*s
rhos_CO2 = CO2_data(:, 3)'.*MW_CO2./1000; %kg/m^3

%interpolate
rho_CO2 = interp1(P_CO2, rhos_CO2, Pcrit, 'spline'); %kg/m^3
mu_CO2 = interp1(P_CO2, mus_CO2, Pcrit, 'spline'); %Pa*s

[vt, ~, ~] = solve_vt(D, mu_CO2, rho_IL, rho_CO2, g); %terminal velocity in m/s


%To solve for critical pressure and x1
%[K, x1crit, Pcrit] = solve_Pcrit(T, R, Tm, dHfus, dH, dS);

%To solve for terminal velocity
%[vt, Re, Cd] = solve_vt(D, visc, rho_IL, rho_CO2, g);

end

function [K, x1crit, Pcrit] = solve_Pcrit(T, R, Tm, dHfus, dH, dS)
%Solves for Pcrit, depends on T

dG = dH - T*dS;
K = exp(-dG/(R*T));

x1crit = exp((-dHfus/R).*((1/Tm) - (1./T)));
Pcrit = (1-x1crit)./(x1crit.*K);
end

function [vt, Re, Cd] = solve_vt(D, visc, rho_IL, rho_CO2, g)
%To solve for terminal velocity
%needs rhos, CO2 viscosity, and diameter

vts = 0:0.01:20; %m/s
Re_vals = Reynolds(vts, D, visc); %vt is terminal velocity
as_mat = dragconstants(Re_vals);
Cds = as_mat(:,1)' + as_mat(:,2)'./Re_vals + as_mat(:,3)'./(Re_vals.^2);

%minimize to find vt
func = 4*rho_IL*g - (3*Cds.*rho_CO2.*(vts.^2));
minfunc = abs(func);
index = (minfunc == min(minfunc));

%final values
vt = vts(index);
Re = Re_vals(index);
Cd = Cds(index);
end

function [Re] = Reynolds(u, D, visc)
%inputs are velocity, diameter, and velocity

Re = u.*D./visc;
end

function [as_mat] = dragconstants(Re_vals)
%drag coefficient equation constants, dependent on Re

as_mat = zeros(length(Re_vals), 3);
for i = 1:length(Re_vals)
    Re = Re_vals(i);
    if Re < 0.1
        as = [0 24 0];
    elseif Re < 1
        as = [3.69 22.73 0.0903];
    elseif Re < 10
        as = [1.222 29.1667 -3.8889];
    elseif Re < 100
        as = [0.6167 46.6 -116.67];
    elseif Re < 1000
        as = [0.3644 98.33 -2778];
    elseif Re < 5000
        as = [0.357 148.62 -47500];
    elseif Re < 10000
        as = [0.46 -490546 578700];
    elseif Re >= 10000
        as = [0.5191 -1662.5 5416700];
    end
    as_mat(i, (1:3)) = as;
end
end

function [T] = solve_T_from_xcrit(x1crit, dHfus, R, Tm)

T = 1./(-(log(x1crit)./(-dHfus/R)) + (1/Tm));
end

function [sol] = co2_solubility(T, P)
%Interpolate experimental data to get CO2 solubility, T in C, P in bar

%%Solubility, mol CO2/IL
P60 = [0.032, 0.077, 0.078, 0.102, 0.202, 0.424, 0.708, 0.951];
sol60 = [0.13, 0.29, 0.5, 0.68, 0.87, 0.91, 0.95, 0.97];

P70 = [0.059, 0.099, 0.119, 0.158, 0.336, 0.699, 0.994];
sol70 = [0.08, 0.27, 0.56, 0.82, 0.86, 0.9, 0.91];

P80 = [0.104, 0.145, 0.201, 0.293, 0.553, 0.787, 0.933];
sol80 = [0.11, 0.22, 0.43, 0.67, 0.76, 0.77, 0.78];

if T == 60
    sol = interp1(P60, sol60, P, 'spline');
elseif T == 70
    sol = interp1(P70, sol70, P, 'spline');      
elseif T == 80
    sol = interp1(P80, sol80, P, 'spline');
else
    T60 = interp1(P60, sol60, P, 'spline');
    T70 = interp1(P70, sol70, P, 'spline');
    T80 = interp1(P80, sol80, P, 'spline');
    sol = interp1([60, 70, 80], [T60, T70, T80], T, 'linear', 'extrap');
end

end
 
function [rho] = IL_density(T)
 
T_range = [10, 20, 22, 30, 40, 50, 60, 70];
density = [1.154, 1.147, 1.145, 1.140, 1.134, 1.128, 1.121, 1.115];
rho = interp1(T_range, density, T, 'spline');
end

