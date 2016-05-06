function CHBE402_Project_Pressure_final2
%constants
g = 9.81; %m/s^2
R = 8.314; %J/molK
%Rbar = 8.3144598e-2; % m^3 * bar / K*kmol
Tm = 166 + 273.15; %K, MP of IL
dHfus = -19900; %J/mol
dH = -58500; %J/mol
dS = -143; %J/mol
MW_CO2 = 44.01; %g/mol
MW_IL = 264.39; %g/mol
MW_Comp = MW_CO2 + MW_IL; %g/mol
D = 0.002; %m
%Dab = .1e-4; %m^2/s at 1bar
Ad = 4*pi*(D/2)^2;
V_drop_eff = (4/3)*pi*(D)^3; %m^3, effective vol of spaced droplet

%Obtain CO2 data and set up interpolation funcs
[Ps_CO2, mus_CO2, rhos_CO2] = CO2_data;
CO2_densityfun = @(P)(interp1(Ps_CO2, rhos_CO2, P, 'spline')); %kg/m^3
CO2_viscosityfun = @(P)(interp1(Ps_CO2, mus_CO2, P, 'spline')); %Pa*s


% Set Initial Conditions
zn = 0.35; % Molefraction of carbon containing components
T = 434; %K
Ac = 5; % m^2
deltaH = 0.005; %m

sphere_packing = 0.52;

N_drops = deltaH*Ac*sphere_packing/V_drop_eff; % Check after pressure gradient is found
rho_IL = IL_density(T);

% Overall Mass Balance
Vapor_out_mass = 596000 / 3600; %kg/s
Vapor_out_molar = Vapor_out_mass / (MW_CO2); % kmol/s
Carbon_in_molar = Vapor_out_molar / (1-zn); %kmol/s
Total_in_molar = 2*Carbon_in_molar;
disp(Total_in_molar*MW_Comp);

ni = Carbon_in_molar; %kmol/s, Equimolar IL and Carbon (ni is an IL molar flow rate)

% Start of Flow chart
% Bottom-most stage
%%%%%%%%
[K, ~, Pcrit] = solve_Pcrit_and_K(T, R, Tm, dHfus, dH, dS);
Pcrit = 1.37;
Ps = zn / (K*(1-zn)); % Surface pressure at bottom of tower (in bar)

rhoS = CO2_densityfun(Ps); %kmol/m^3
rhoN = CO2_densityfun(Pcrit); %kmol/m^3
rhoN_mass = rhoN*MW_CO2; %kg/m^3

viscosity = CO2_viscosityfun(Pcrit);
vt = solve_vt(D, viscosity, rho_IL, rhoN_mass, g);
Re = Reynolds(vt, D, viscosity, rhoN_mass);

Dab = get_Dab(Pcrit);
hm = hmfunc(Dab, Re, viscosity, rhoN_mass, deltaH);

% Solve for zn-1 from convective mass transfer
%zMinus1 = zn + (hm*Ad*N_drops*(rhoS-rhoN))/ni
rhoN = hm*Ad*N_drops*rhoS/(Ac*deltaH + hm*Ad*N_drops);
zMinus1 = zn + Ac*deltaH*rhoN/ni;
%deltaH = ni*(zMinus1 - zn)/(Ac*rhoN)


% Solve for vapor molar flow rate at bottom
%Vn = ni * (zMinus1 - zn); % mol/s
%Q = ni * (zMinus1 - zn)/rhoN;

disp('    Ps      rhoS    hm          rhoN    zMinus1  P      viscosity   vt      Re')
fprintf('\t%.3f \t%.3f \t%.6f \t%.3f \t%.3f \t%.3f \t%.5f \t%.2f \t%.f\n', Ps, rhoS, hm, rhoN, zMinus1, Pcrit, viscosity, vt, Re)
disp(' ')
%disp(' ')

% End of bottom stage
disp('    Ps      rhoS    hm          rhoN    zMinus1  P      viscosity   vt      Re')
% While loop for the next stages 
counter = 0;
while zMinus1 < .99
   counter = counter + 1;
   
   % Record values brought from previous iteration
   z = zMinus1;
   rhoMinus1 = rhoN; %kmol/m^3
   rhoMinus1_mass = rhoMinus1*MW_CO2; %kg/m^3
   
   Ps = z/(K*(1-z));
   rhoS = CO2_densityfun(Ps); %kmol/m^3
   
   hm = hmfunc(Dab, Re, viscosity, rhoMinus1_mass, deltaH);
   
   %Helen's new code:
   rhoN = (hm*Ad*sphere_packing*rhoS/V_drop_eff)/(1+ (hm*Ad*sphere_packing/V_drop_eff));
   zMinus1 = z + hm*Ad*N_drops*(rhoS-rhoN)/ni;
   
   % Old code: 
   % Use convection equation to solve for the density at infinity
   %rhoN = rhoS / (1+Q/(hm*Ad*N_drops)); %kmol/m^3
   % Use the new density to solve for the zMinus1
   %zMinus1 = (Q*rhoN)/(ni*MW_CO2) + z;
   
   deltaRhoVec(counter) = rhoMinus1 - rhoN;
   hmvec(counter) = hm;
   
   rhoN_mass = rhoN*MW_CO2; %kg/m^3

   % Solve for P
   zerofun = @(P)(CO2_densityfun(P) - rhoN);
   P = fzero(zerofun,1);
   Pvec(counter) = P;
   
   % Re-evaluate viscosity
   viscosity = CO2_viscosityfun(P);
   vt = solve_vt(D, viscosity, rho_IL, rhoN_mass, g);
   Re = Reynolds(vt,D,viscosity,rhoN_mass);
   Dab = get_Dab(P);
   %Q = ni * (zMinus1 - zn) / rhoN;
   
   fprintf('%.f \t%.3f \t%.3f \t%.6f \t%.3f \t%.3f \t%.3f \t%.5f \t%.2f \t%.f\n', counter, Ps, rhoS, hm, rhoN, zMinus1, P, viscosity, vt, Re)
end

disp('       Ps      rhoS    hm          rhoN    zMinus1  P      viscosity   vt      Re')

%deltaT = ni * (zn - z) / ((Ad*N_drops) * sum(hmvec.*deltaRhoVec));
%height = deltaT*(counter-1)*vt
height = counter*deltaH;
N_drops

figure(1)
plot(linspace(0,height,counter-5),Pvec(1:end-5),'linesmoothing','on','linewidth',1.5);
xlabel('Tower Height (m)')
ylabel('Pressure (bar)')
%cut off last iteration because z > 1
%Hvec = Hvec(1:end-1);
%Pvec = Pvec(1:end-1);

%plot(Hvec, Pvec)
%xlabel('Height (m)')
%ylabel('Pressure (bar)')

end

function [out] = hmfunc(Dab,Re,viscosity,rho,L)
    % L is the deltaH at a stage
    % V is terminal velocity in most cases
    nu = viscosity / rho;
    Sc = nu/Dab;
    Sh = 2 + 0.6 * sqrt(Re) * Sc^(1/3);
    out = Dab*Sh/L;
end

function [Ps_CO2, mus_CO2, rhos_CO2] = CO2_data

%Open file for CO2 data: viscosity and density
%Based on T = 415.27
CO2_fileID = fopen('CO2_data.txt');
CO2_data = fscanf(CO2_fileID, '%f', [13, 58])';

Ps_CO2 = CO2_data(:, 2)'; %bar
mus_CO2 = CO2_data(:, 12)'; %Pa*s
rhos_CO2 = CO2_data(:, 3)'./1000; %kmol/m^3

end

function [K, x1crit, Pcrit] = solve_Pcrit_and_K(T, R, Tm, dHfus, dH, dS)
%Solves for Pcrit, depends on T

dG = dH - T*dS;
K = exp(-dG/(R*T));

x1crit = exp((-dHfus/R).*((1/Tm) - (1./T)));
Pcrit = (1-x1crit)./(x1crit.*K);
end

function [vt, Re, Cd] = solve_vt(D, visc, rho_IL, rho_CO2m, g)
%To solve for terminal velocity
%needs rhos, CO2 viscosity, and diameter

vts = 0:0.01:50; %m/s
Re_vals = Reynolds(vts, D, visc, rho_CO2m); %vt is terminal velocity
as_mat = dragconstants(Re_vals);
Cds = as_mat(:,1)' + as_mat(:,2)'./Re_vals + as_mat(:,3)'./(Re_vals.^2);

%minimize to find vt
func = 4*rho_IL*g*D - (3*Cds.*rho_CO2m.*(vts.^2));
minfunc = abs(func);
index = (minfunc == min(minfunc));

%final values
vt = vts(index);
Re = Re_vals(index);
Cd = Cds(index);
end

function out = get_Dab(P)
    out = .2e-4 * 1/(P/.506);
end

function [Re] = Reynolds(u, D, visc, rho)
%inputs are velocity, diameter, and velocity

Re = u.*rho.*D./visc;
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

function [rho] = IL_density(T)
% T in K

T = T - 273.15; %convert to C
T_range = [10, 20, 22, 30, 40, 50, 60, 70];
density = [1.154, 1.147, 1.145, 1.140, 1.134, 1.128, 1.121, 1.115];
rho = interp1(T_range, density, T, 'spline'); % g/cm^3
rho = rho*1000; %kg/m^3
end