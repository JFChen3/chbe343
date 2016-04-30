function CHBE402_Project_Helen

%constants
g = 9.81; %m/s^2
R = 8.314; %J/molK
Tm = 166 + 273.15; %MP of IL
dHfus = -19900; %J/mol
dH = -58500; %J/mol
dS = -143; %J/mol*K

%Test calculation
T = 423;

%To solvce for critical pressure and x1
[K, x1crit, Pcrit] = solve_Pcrit(T, R, Tm, dHfus, dH, dS);

Pcrit

%To solve for terminal velocity
[vt, Re, Cd] = solve_vt(D, visc, rho_il, rho_co2, g);

end

function [K, x1crit, Pcrit] = solve_Pcrit(T, R, Tm, dHfus, dH, dS)
%Solves for Pcrit, depends on T

dG = dH - T*dS;
K = exp(-dG/(R*T));

x1crit = exp((-dHfus/R).*((1/Tm) - (1./T)));
Pcrit = (1-x1crit)./(x1crit.*K);
end

function [vt, Re, Cd] = solve_vt(D, visc, rho_il, rho_co2, g)
%To solve for terminal velocity
%needs rhos, CO2 viscosity, and diameter

vts = 0:0.01:20; %m/s
Re_vals = Reynolds(vts, D, visc); %vt is terminal velocity
as_mat = dragconstants(Re_vals);
Cds = as_mat(:,1)' + as_mat(:,2)'./Re_vals + as_mat(:,3)'./(Re_vals.^2);

%minimize to find vt
func = 4*rho_il*g - (3*Cds.*rho_co2.*(vts.^2));
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