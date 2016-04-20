function air_conditioner

C_power = 3.6; %kWh/rev
B_meter = 14.4; %kWh/rev
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Low flow rate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = [72.6, 73.3, 66.4, 61.9, 79.2, 55.7, 47.6, 50.6, 52.5]; %F
T_use = convert_to_Kelvin(T);
H_air = [67.45, 44.19]; %Enthalpy of air going into and out of evaporator, kJ/kg-dry
H_condensate = 37.80; %kJ/kg, assumed to be saturated liquid at exit wet bulb
H_refrig = [20.045, 35.771]; %Enthalpy of refrigerant going into and out of evaporator, kJ/mol
mass_ratio = [0.013, 0.0065]; %Mass water/mass dry air, inlet and outlet
comp_time = 8.89; %s
blow_time = 231.5; %s
power = calc_power([C_power, B_meter],[comp_time,blow_time]); %kW
ref_flow = 758; %Refrigerant pulse

disp('--------------------------LOW FLOW RATE--------------------------')
%Coefficient of performance
[CoP] = coeff_performance(ref_flow, H_refrig, power);
fprintf('Coefficient of performance: %f\n',CoP)

anemometer = 1.22;
humidity_in = 0.638;
humidity_out = 0.973;
pressure_1 = 65 + P_atm;
pressure_2 = 195 + P_atm;
condense_flow = 362/378; %mL/s

calc_all(H_air, H_condensate, mass_ratio, anemometer, condense_flow, T_use(1), ref_flow, H_refrig, [C_power, B_meter], [comp_time,blow_time])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%High flow rate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = [72.8, 73.3, 66.3, 63.9, 80.8, 61.3, 56.4, 57.8, 59.8]; %F
T_use = convert_to_Kelvin(T);
H_air = [67.78, 56.99]; %kJ/kg-dry
H_condensate = 37.80; %kJ/kg, assumed to be saturated liquid at exit wet bulb
H_refrig = [20.129, 35.879]; %kJ/mol
mass_ratio = [0.013, 0.010]; %Mass water/mass dry air, inlet and outlet
comp_time = 8.53; %s
blow_time = 144.5; %s
power = calc_power([C_power, B_meter],[comp_time,blow_time]); %kW
ref_flow = 843; %Refrigerant pulse

disp('--------------------------HIGH FLOW RATE--------------------------')
%Coefficient of performance
[CoP] = coeff_performance(ref_flow, H_refrig, power);
fprintf('Coefficient of performance, low flow rate: %f\n',CoP)

anemometer = 1.44;
humidity_in = 0.650;
humidity_out = 0.955;
pressure_1 = 73 + P_atm;
pressure_2 = 203 + P_atm;
condense_flow = 229/270; %mL/s

calc_all(H_air, H_condensate, mass_ratio, anemometer, condense_flow, T_use(1), ref_flow, H_refrig, [C_power, B_meter], [comp_time,blow_time])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%Expansion valve data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P1 = 50
P2_50 = [55,50:-1:32]; %psig
flow_rate_50 = [0,0,0.8,1.2,1.7,2.2,2.6,3.0,3.2,3.5,3.8,3.9,4.0,4.2,4.3,4.3,4.4,4.4,4.4,4.4]; %scfm

%P1 = 15
P2_15 = 50:-1:34;
flow_rate_15 = [0,0,0,0.8,0.9,1.0,1.2,1.2,1.4,1.4,1.5,1.6,1.6,1.6,1.7,1.7,1.7];

plot(P2_50,flow_rate_50,'.-',P2_15,flow_rate_15,'o-')
xlabel('Flow Rate (scfm)')
ylabel('P_2 (psig)')
legend('P_1 = 50 psig','P_2 = 15 psig')

end

function [T_K] = convert_to_Kelvin(T_F)
%Convert Fahrenheit to Kelvin

T_K = (T_F+32).*(5/9) + 273;

end

function calc_all(H_air, H_condensate, mass_ratio, anemometer, condense_flow, T_in, refrig_pulse, H_refrig, power_per_rev, rev_time)

%Mass flow rates of air and condensate
[~, m_in_dry, ~, m_out_dry, m_water] = mass_flow_air(anemometer, mass_ratio, condense_flow, T_in);

%Enthalpy of air streams
H_in = m_in_dry*H_air(1); %kJ/s
H_out = m_out_dry*H_air(2); %kJ/s

%Energy balance to find heat loss
refrig_flow = molar_flow_rate_refrig(refrig_pulse); %mol/s
Q_ref = refrig_flow*(H_refrig(2) - H_refrig(1)); %kJ/s
Q_air = H_out - H_in;
Q_condensate = m_water*H_condensate;

Q_loss = abs(Q_condensate + Q_air + Q_ref);

%Print results
fprintf('Specific enthalpy, inlet air: %.2f kJ/kg-dry\n',H_air(1))
fprintf('Total enthalpy, inlet air: %.2f kJ/s\n',H_in)
fprintf('Specific enthalpy, outlet air: %.2f kJ/kg-dry\n',H_air(2))
fprintf('Total enthalpy, outlet air: %.2f kJ/s\n',H_out)
fprintf('Specific enthalpy, condensate: %.2f kJ/kg\n',H_condensate)
fprintf('Total enthalpy, condensate: %.4f kJ/s\n',Q_condensate)
fprintf('Heat loss: %.2f kJ/s\n',Q_loss)

end

function [CoP] = coeff_performance(refrig_pulse, H_refrig, power)

refrig_flow = molar_flow_rate_refrig(refrig_pulse); %mol/s
Qc = refrig_flow*(H_refrig(2) - H_refrig(1)); %kJ/s
CoP = Qc/power; %dimensionless

end

function [eff] = compressor_efficiency(W, T1, P1, P2)
%Compressor efficiency

Ws = Cp*T1*((P2/P1)^(R/Cp) - 1); %Sourced from thermo textbook
eff = W/Ws;

end

function [m_in, m_in_dry, m_out, m_out_dry, m_water] = mass_flow_air(anemometer, mass_ratio, condense_flow, T_in)
%Mass flow rate of air into evaporator, air and water out of evaporator

mass_ratio_in = mass_ratio(1);
mass_ratio_out = mass_ratio(2);

A = 14*22*(0.0254)^2; %m^2, cross sectional area of duct
v = 0.8*anemometer_cal(anemometer)*0.3048; %m/s, average velocity
Q_air = A*v;

%Inlet air stream
density = density_air(mass_ratio_in, T_in);
m_in = density*Q_air;

%Exiting water stream
m_water = condense_flow*0.001; %Convert mL/s to kg/s

%Outlet air stream
m_out = m_in - m_water;

%Dry air mass flow rates
m_in_dry = m_in*(1-mass_ratio_in);
m_out_dry = m_out*(1-mass_ratio_out);

end

function [density] = density_air(mass_ratio, T)
%Calculate air density

P = 101325; %1atm in Pa
R = 8.3144598; %m^3Pa/Kmol
M_air = 0.02897; %Molar mass, kg/mol
M_water = 0.01802;

xw = mass_ratio_to_x(mass_ratio);
density = (P*M_air/(R*T))*(1-xw*(1-(M_water/M_air))); %kg/m^3

end

function xw = mass_ratio_to_x(mass_ratio)
%mass_ratio is mass water/1 unit mass dry air
M_air = 28.97;
M_water = 18.02;

mol_air = 1/M_air;
mol_water = mass_ratio/M_water;
mol_total = mol_air + mol_water;

xw = mol_water/mol_total;

end

function [flow_rate] = molar_flow_rate_refrig(pulse)

flow_rate = (pulse*0.03204*0.000001052*1184)/0.08647; %mol/s
                    % gph  m3/s       *kg/m3=kg/s  /(kg/mol) = mol/s

end

function [power] = calc_power(power_per_rev, rev_time)
%Compute total power input, power_per_rev in kWh/rev, rev_time in s

power = sum(power_per_rev./(rev_time/3600));

end

function [velocity] = anemometer_cal(voltage)
%Centerline velocity from anemometer voltage

velocity = 12.02*voltage^2 - 13.627*voltage + 3.7286;

end