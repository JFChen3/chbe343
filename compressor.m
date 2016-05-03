function compressor

rate = 1.3*(10^6)*0.453592; %kg/hr
Pi = 60*6894.76 + 101325; %Pa
Ti = (85-32)*(5/9) + 273; %C
Tf = Ti;
Pf = 2200*6894.76 + 101325; %Pa
R = 8.314; %m^3*Pa/K*mol
Cp = 846; %J/kgK
Cv = 657; %J/kgK
MW = (12 + 16 + 16)/1000; %kg/mol
ratio = Pf/Pi;
eff = 0.65;

k = Cp/Cv;

[minpower, T2_min] = calcminpower(R, Ti, k, MW, Pf, Pi, rate);
powerUS = minpower*0.00134102; %hp

powerUS

part2(rate, Pi, Ti, Tf, Pf, R, Cp, k, MW, ratio, eff)

end

function [minpower, T2_min] = calcminpower(R, Ti, k, MW, Pf, Pi, rate)

work = (R*Ti*k/(MW*(k-1)))*(((Pf/Pi)^((k-1)/k))-1); %J/kg
T2_min = Ti*((Pf/Pi)^((k-1)/k));
minpower = work*rate/(60*60); %J/s or W
end


function part2(rate, Pi, Ti, Tf, Pf, R, Cp, k, MW, ratio, eff)

ratioper = ratio^(1/3);

%compressor 1
Pf1 = Pi*ratioper;
work1 = (R*Ti*k/(MW*(k-1)))*(((Pf1/Pi)^((k-1)/k))-1); %J/kg
Tf1 = Ti*((Pf1/Pi)^((k-1)/k)) + ((1-eff)/eff)*(work1/Cp);

workper = work1/0.65; %J/kg
workperUS = (workper*0.453592)/1.355818; %ftlb/lb
powerper = workper*rate/(60*60); %J/s
powerperUS = powerper*0.00134102; %hp
powerptotUS = powerperUS*3; %hp

%FOR COOLING
%water props
Tw = (68 - 32)*(5/9) + 273; %K
Cpw = 4182; %J/kgK

waterflowper = rate*Cp*(Tf1-Ti)/(Cpw*(Ti-Tw)); %kg/hr
waterflowperUS = waterflowper/(0.453593*60*60); %lb/s
waterflowtotUS = waterflowperUS*3;

Tf1US = (Tf1 - 273)*(9/5) + 32;

powerper
waterflowperUS
waterflowtotUS
Tf1US
end
