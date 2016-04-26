function CHBE402_Project_Helen

%constants
g = 9.81; %m/s^2

%To solve for terminal velocity
%needs rhos, CO2 viscosity, and diameter
vt = 0:0.1:20; %m/s
Re = Reynolds(vt, D, v); %vt is terminal velocity
as = dragconstants(Re);
Cd = as(1) + as(2)./Re + as(3)./(Re.^2);

%vt = sqrt(4*rho_il*g/(3*Cd*rho_co2);

func = 4*rho_il*g - (3*Cd.*rho_co2.*(vt.^2));
minfunc = abs(func);
index = (minfunc == min(minfunc));
final_vt = vt(index);
final_Re = Re(index);
final_Cd = Cd(index);

end

function [Re] = Reynolds(u, D, v)
%inputs are velocity, diamter, and velocity

Re = u.*D./v;
end

function [as] = dragconstants(Re)
%drag coefficient equation constants, dependent on Re

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

end