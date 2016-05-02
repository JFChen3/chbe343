function CHBE402_Project_Helen_Cleanup
%Current state of affairs: Helen cleans up variable names + units
close all

%constants
g = 9.81; %m/s^2
R = 8.314; %J/molK
Rbar = 8.3144598e-2; % m^3 * bar / K*kmol
Tm = 166 + 273.15; %K, MP of IL
dHfus = -19900; %J/mol
dH = -58500; %J/mol
dS = -143; %J/mol
MW_CO2 = 44.01; %g/mol
MW_IL = 264.39; %g/mol
MW_Comp = MW_CO2 + MW_IL; %g/mol
D = 0.002; %m
Dab = .1e-4; %m^2/s at 1bar

%Obtain CO2 data and set up interpolation funcs
[Ps_CO2, mus_CO2, rhos_CO2, vols_CO2] = CO2_data(MW_CO2);
CO2_densityfun = @(P)(interp1(Ps_CO2, rhos_CO2, P, 'spline')); %kg/m^3
CO2_viscosityfun = @(P)(interp1(Ps_CO2, mus_CO2, P, 'spline')); %Pa*s

%Bottom stage: set xi, find T, Pcrit, K, sol
xi = 0.75; % Mole fraction of IL
T = 273.15 + 125; %K
error = 100;
tol = 0.001;
dT = 0.01; %K, temp step
while error > tol
    T = T + dT;
    [K, ~, Pcrit] = solve_Pcrit_and_K(T, R, Tm, dHfus, dH, dS); %Pcrit in bar
    sol = CO2_solubility((T), Pcrit);
    z = solve_z(K, Pcrit); %NOTE: Pcrit converted to Pa?!!!!!!!!?
    error = abs(sol-z)/sol;
end

fprintf('Temperature: %.2f k\n', T);

%Get CO2 properties at the T we found
rho_CO2 = CO2_densityfun(Pcrit); %kg/m^3
mu_CO2 = CO2_viscosityfun(Pcrit); %Pa*s
rho_IL = IL_density(T);
[vt, ~, ~] = solve_vt(D, mu_CO2, rho_IL, rho_CO2, g); %terminal velocity in m/s



%%%%%%%%%%%%%% START CODE THAT NEEDS UPDATING TO NEW K DEF %%%%%%%%%%%%%%%%%

% Molar Flow out of bottom of column
% Overall Mass Balance
Vapor_out_mass = 596000 / 3600; %kg/s
Vapor_out_molar = Vapor_out_mass / (MW_CO2); % kmol/s
Feed_CO2_molar = Vapor_out_molar * (2-xi); %kmol/s %%%%%%%%%%% This looks suspicious pls check
Feed_ionic_molar = Feed_CO2_molar;
Feed_mass = Feed_CO2_molar * (MW_CO2) + Feed_ionic_molar * (MW_IL); %kg/s
Bottoms_out_molar = Feed_ionic_molar + Feed_CO2_molar * (1-xi);
Bottoms_out_mass = Bottoms_out_molar * (xi*MW_IL + (1-xi)*MW_CO2); %kg/s

Liquid_Molar = Bottoms_out_molar; % kmol/s

% Choose a height increment
Area = 10; %m^2
total_drops = Area / (2*D);
dH = 2; %m
dT = dH/vt; % If v = m/s, deltaT = s

% Surface Area of the Drop (consider constant)
A_drop = 4 * pi * (D/2)^2; % M^2

% Stage 1, Helen updated mol fracs
xi_bot = xi;
xcomp_bot = K*xi_bot*Pcrit;
xco2_bot = sol - xcomp_bot;
% xco2_bot = xco2_from_xi(xi_bot,K);
% xcomp_bot = 1 - xi_bot - xco2_bot;
%%%%%%%%%%%%%%%%%%%%%%

% Functions needed
Re = Reynolds(vt,D,CO2_viscosityfun(Pcrit),CO2_densityfun(Pcrit));
hm = hmfunc(Dab,Re,CO2_viscosityfun(Pcrit),CO2_densityfun(Pcrit),dH);

%currently, rho_as is density of CO2 dissolved in liquid
% and rho_inf is density of CO2 at that stage
rho_CO2s = rho_IL / (xi_bot*MW_IL + xcomp_bot*MW_Comp + xco2_bot*MW_CO2) * xco2_bot * MW_CO2;
Vapor_flow_n = total_drops*hm*A_drop*dT*(rho_CO2s - CO2_densityfun(Pcrit));

% Flow rate of Ln-1 * xn-1
Lcomp = Liquid_Molar*(xcomp_bot+xco2_bot) - Vapor_flow_n;

% Get xCO2 + xComplex of the stage above (both carbon containing species);
x_combined = Lcomp / Liquid_Molar;

% xi = 1 - x_combined;
% xco2 = xco2_from_xi(xi,K);
% xcomp = 1 - xi - xco2;
xco2 = xco2_bot;
xcomp = xcomp_bot;

% Create the loop to iterate up the column
counter = 0;
Pvec = zeros(80);
Hvec = zeros(80);
while counter < 80;
    counter = counter + 1;
    
    % VnMinus is coming up from the layer below
    if counter == 1;
        VnMinus = Vapor_flow_n;
        Hvec(counter) = dH;
    else
        VnMinus = Vn;
        Hvec(counter) = Hvec(counter-1) + dH;
    end
    
    % Find the pressure of the next stage using the concentration of IL
    % Using density from stage below
    Sol = xcomp + xco2;
    
    % Numerically solve to find the P which delivers the necessary solubility
    zerofun = @(P)(CO2_solubility(T, P) - Sol);
    P = fzero(zerofun, 0.5); %bar   
    Pvec(counter) = P; %Record pressure at this stage
    
    %Update substance properties
    rho_IL = IL_density(T); %kg/m^3
    rho_CO2 = CO2_densityfun(P); % kg/m^3
    mu_CO2 = CO2_viscosityfun(P); %Pa*s
    
    Re = Reynolds(vt, D, mu_CO2, rho_CO2);
    hm = hmfunc(Dab, Re, mu_CO2, rho_CO2, dH);
        
    % Solve for vapor flow rate out of the stage
    rho_CO2s = (rho_IL / (xi*MW_IL + xcomp*MW_Comp + xco2*MW_CO2)) * xco2 * MW_CO2;
    
    Vn = total_drops*hm*A_drop*dT*(rho_CO2s - CO2_densityfun(P));
    Vn = Vn * 30000;
    deltaVn = Vn - VnMinus;
    
    % Solve for the Liquid flow rate coming into the stage
    LcompMinus = Liquid_Molar*(xcomp + xco2) + deltaVn;
       
    % Get new combined molefraction of carbon containing species
    x_combined = LcompMinus / Liquid_Molar;
    
    xi = 1 - x_combined;
    xco2 = xco2_from_xi(xi,K);
    xcomp = 1 - xco2 - xi;
end

plot(Hvec,Pvec)
xlabel('Height (m)')
ylabel('Pressure (bar)')

end

function [out] = xco2_from_xi(xi,K)
numer = 1 - xi;
denom = (xi*K) * (1 + 1/(xi*K));
out = numer / denom;
end

function [K, x1crit, Pcrit] = solve_Pcrit_and_K(T, R, Tm, dHfus, dH, dS)
%Solves for Pcrit, depends on T

dG = dH - T*dS;
K = exp(-dG/(R*T));

x1crit = exp((-dHfus/R).*((1/Tm) - (1./T)));
Pcrit = (1-x1crit)./(x1crit.*K);
end

function [vt, Re, Cd] = solve_vt(D, visc, rho_IL, rho_CO2, g)
%To solve for terminal velocity
%needs rhos, CO2 viscosity, and diameter

vts = 0:0.01:50; %m/s
Re_vals = Reynolds(vts, D, visc, rho_CO2); %vt is terminal velocity
as_mat = dragconstants(Re_vals);
Cds = as_mat(:,1)' + as_mat(:,2)'./Re_vals + as_mat(:,3)'./(Re_vals.^2);

%minimize to find vt
func = 4*rho_IL*g*D - (3*Cds.*rho_CO2.*(vts.^2));
minfunc = abs(func);
index = (minfunc == min(minfunc));

%final values
vt = vts(index);
Re = Re_vals(index);
Cd = Cds(index);
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

function [sol] = CO2_solubility(T, P)
%Interpolate experimental data to get CO2 solubility, T in K, P in bar

T = T - 273.15; %convert to C

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
 
function [sol] = solve_z(K, P)
%Equilibrium equation, returns (xc + xco2)

sol = K.*P./(1 + K.*P);
end

function [rho] = IL_density(T)
% T in K

T = T - 273.15; %convert to C
T_range = [10, 20, 22, 30, 40, 50, 60, 70];
density = [1.154, 1.147, 1.145, 1.140, 1.134, 1.128, 1.121, 1.115];
rho = interp1(T_range, density, T, 'spline'); % g/cm^3
rho = rho*1000; %kg/m^3
end

function [Ps_CO2, mus_CO2, rhos_CO2, vols_CO2] = CO2_data(MW_CO2)

%Open file for CO2 data: viscosity and density
%Based on T = 415.27
CO2_fileID = fopen('CO2_data.txt');
CO2_data = fscanf(CO2_fileID, '%f', [13, 193])';

Ps_CO2 = CO2_data(:, 2)'; %bar
mus_CO2 = CO2_data(:, 12)'; %Pa*s
rhos_CO2 = CO2_data(:, 3)'.*MW_CO2./1000; %kg/m^3
vols_CO2 = CO2_data(:,4); % Molar Volume (m^3/mol)

end

function [out] = hmfunc(Dab,Re,viscosity,rho,L)
    % L is the deltaH at a stage
    % V is terminal velocity in most cases
    nu = viscosity / rho;
    Sc = nu/Dab;
    Sh = 2 + 0.6 * sqrt(Re) * Sc^(1/3);
    out = Dab*Sh/L;
end
