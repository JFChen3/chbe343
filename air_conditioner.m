function air_conditioner

C_power = 3.6; %kWh/rev
B_meter = 14.4; %kWh/rev
pulse_scale = 0.03204; %Convert refrigerant flow to gph (who the fuck uses gallons)
P_atm = 14.7;

% TC 1 – Intake air temperature
% TC 2 – Intake wet bulb temperature
% TC 3 – Intake dry bulb temperature
% TC 4 – Condensate exiting evaporator
% TC 5 – High side refrigerant (liquid entering expansion valve from compressor)
% TC 6 – Low side refrigerant (vapor leaving the evaporator coils & returning to compressor)
% TC 7 – Exit wet bulb temperature (between evaporator and blower)
% TC 8 – Exit dry bulb temperature (between evaporator and blower)
% TC 9 – Exit temperature (at the blower)

%%Low flow rate data
T = [72.6, 73.3, 66.4, 61.9, 79.2, 55.7, 47.6, 50.6, 52.5]; %F
H_air = [20.5, 19.7]; %Enthalpy of air going into and out of evaporator, Btu/lb-dry
H_refrig = [20.045, 35.771]; %Enthalpy of refrigerant going into and out of evaporator, kJ/mol
comp_time = 8.89; %s
blow_time = 231.5; %s
power = calc_power([C_power, B_meter],[comp_time,blow_time]); %kW
ref_flow = 758; %Refrigerant pulse

%Coefficient of performance
[CoP] = coeff_performance(ref_flow, H_refrig, power);
fprintf('Coefficient of performance, high flow rate: %f\n',CoP)

Anemometer = 1.22;
humidity_in = 0.638;
humidity_out = 0.973;
pressure_1 = 65 + P_atm;
pressure_2 = 195 + P_atm;
condense_flow = 362/378; %mL/s

%%High flow rate data
T = [72.8, 73.3, 66.3, 63.9, 80.8, 61.3, 56.4, 57.8, 59.8]; %F
H_air = [20.5, 21.5];
H_refrig = [20.129, 35.879];
comp_time = 8.53; %s
blow_time = 144.5; %s
power = calc_power([C_power, B_meter],[comp_time,blow_time]); %kW
ref_flow = 843; %Refrigerant pulse

%Coefficient of performance
[CoP] = coeff_performance(ref_flow, H_refrig, power);
fprintf('Coefficient of performance, low flow rate: %f\n',CoP)

Anemometer = 1.44;
humidity_in = 0.650;
humidity_out = 0.955;
pressure_1 = 73 + P_atm;
pressure_2 = 203 + P_atm;
condense_flow = 229/270;

%%Expansion valve data
%P1 = 50
P2_50 = [55,50:-1:32]; %psig
P2_50a = P2_50 + 14.7; % psia
flow_rate_50 = [0,0,0.8,1.2,1.7,2.2,2.6,3.0,3.2,3.5,3.8,3.9,4.0,4.2,4.3,4.3,4.4,4.4,4.4,4.4]; %scfm

%P1 = 15
P2_15 = 50:-1:34;
P2_15a = P2_15 + 14.7; % psia
flow_rate_15 = [0,0,0,0.8,0.9,1.0,1.2,1.2,1.4,1.4,1.5,1.6,1.6,1.6,1.7,1.7,1.7];

psat = r22_sat_pressure(273.15);
title('Flow Rate versus Pressure Differential')
fa_50 = psat - P2_50a;
fa_15 = psat - P2_15a;
prop_fa_50 = fa_50 / psat;
prop_fa_15 = fa_15 / psat;
temp_50 = 273.15 * prop_fa_50;
temp_15 = 273.15 * prop_fa_15;
xlabel('Flow Rate (SCFM)')
ylabel('Psat - P2')

figure(1)
hold on
plot(flow_rate_50,fa_50)
plot(flow_rate_15,fa_15)

figure(2)
plot(P2_50,flow_rate_50,'.-',P2_15,flow_rate_15,'o-')
ylabel('Flow Rate (scfm)')
xlabel('P_2 (psig)')
legend('P_1 = 50 psig','P_1 = 15 psig')

end

function [CoP] = coeff_performance(refrig_pulse, H_refrig, power)

refrig_flow = molar_flow_rate_refrig(refrig_pulse); %mol/s
Qc = refrig_flow*(H_refrig(2) - H_refrig(1)); %kJ/s
CoP = Qc/power; %dimensionless

end

function [flow_rate] = molar_flow_rate_refrig(pulse)

flow_rate = (pulse*0.03204*0.000001052*1184)/0.08647; %mol/s

end

function [power] = calc_power(power_per_rev, rev_time)
%Compute total power input, power_per_rev in kWh/rev, rev_time in s

power = sum(power_per_rev./(rev_time./3600));

end

function ac_performance()

end

function out = r22_sat_pressure(T)
% Antoine equation (temp in K)
% Constants found on Nist Webbook based on Stull 1947
A = 4.36567;
B = 947.577;
C = -14.964;
logP = A - (B/(T+C));
out = 10^(logP); % bar
out = out * 14.5038; % bar to psi

end

function [velocity] = anemometer_cal(pulse)

velocity = 12.02*pulse^2 - 13.627*pulse + 3.7286;

end

velocity = 12.02*pulse^2 - 13.627*pulse + 3.7286;

end
