function chbe402project_interp

%Convenience functions, nothing to see here

%Test calculation
co2_solubility(424-273, 0.4)
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
