function AC_scaleup
% CHANGE THIS TO ARGUMENT
COP = 4;
efficiency = .53;

MW_R22 = 86.47; % g/mol

% Airflow
people = 6; % Estimated people occupying home
area = 5700; % Sq. ft.
af_req = 7.5 * people + .01 * area; % AHRE 62.2 Required airflow in cfm
%disp(af_req)
af_req = af_req * .4719474499999; % Convert to L/s
% Questionable assumption
af_req = af_req * 2.5; 
humidity = .89;

air_temp = 294.261; % 70 degF
air_pressure = 14.7; %psia

P_water_sat = sat_pressure_water(air_temp);
P_water = P_water_sat * humidity;

x_water = P_water / air_pressure;

molar_mass = 28.02 * (1 - x_water) + 18.01528 * x_water; %g/mol

molar_airflow = af_req * 101.33 / (8.314 * 21.1111);
mass_airflow = molar_airflow * molar_mass; % gram / s;
mass_airflow = mass_airflow / 1000; % kg/s

cp_air = 1.009 * (1 - x_water) + 4.18 * x_water; % kJ/kg*K

    function out = sat_pressure_water(T)
        % Antoine equation to find sat. press. of water
        % Coefficients given by Bridgeman and Aldrich
        
        % T in kelvin
        if T >= 273 && T < 304
            A = 5.40221;
            B = 1838.675;
            C = -31.737;
            
        elseif T >= 304 && T < 334
            A = 5.20389;
            B = 1733.926;
            C = -39.485;
            
        elseif T >= 334 && T < 344
            A = 5.0768;
            B = 1659.793;
            C = -45.845;
            
        elseif T >= 344 && T <= 373
            A = 5.08364;
            B = 1663.125;
            C = -45.622;
        end
        
        logP = A - (B/(T+C));
        out = 10^logP; % Pressure in bar
        out = out * 14.5038; % Convert to psi
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

% Given by CoolCalc
Qc = 69500; %BTU/hr

% To change to kJ/s (kW)
Qc = Qc * 1055.06; %J/hr
Qc = Qc / 60; %J/min
Qc = Qc / 60; %J/s (W)
Qc = Qc / 1000; %kJ/s (kW)

% Find the temperature of outgoing air
t_in = air_temp;
t_out = t_in - Qc / (mass_airflow * cp_air);

t_r22_in = 308.706; %K (96 degF)
p_high = r22_sat_pressure(319.817); %psia (116 degF)

t_r22_out = t_out + 2.8; %5 degF above the outlet air temperature
p_low = r22_sat_pressure(t_r22_out - 10.11); %psia (10.11 degC is 18.2 degF)

disp(p_low)

% Enthalpies
% R22 in: P = 262 psia, T = 308.7 K
H_in_22 = 21.074; %kJ/mol

% R22 out: P = 74.9 psia, T = 287.2849 K
H_out_22 = 35.807; %kJ/mol

H_comp_22 = 38.786;

% Choose R-22 enthalpy change in the high flow condition
    function COP_calc(H_in,H_out,MW,H_comp)
        dH = H_out - H_in; %kJ/mol
        
        disp(dH)

        molar_flow = Qc / dH; % mol / s

        mass_flow = molar_flow * MW; % g / s
        
        if exist('H_comp')
            dH_comp = H_comp - H_out;
            Ws = dH_comp * molar_flow;
            W = Ws / efficiency;
            fprintf('Power for Compressor = %.3f kW\n',W)
        end

        compressor_power = Qc/COP; % kJ/s (kW)

        fprintf('Compressor Power: %.3f kW\n',compressor_power)
        fprintf('Refrigerant Flow: %.3f g/s\n',mass_flow)
    end

COP_calc(H_in_22,H_out_22,MW_R22,H_comp_22)

% Part Two, changing the Refrigerant
% 

    function out = prop_sat_pressure(T)
        % T in k
        % P in bar
        if T >= 230.6 && T < 320.7
            A = 3.98292;
            B = 819.296;
            C = -24.417;
        end
        
        logP = A - (B/(T+C));
        P = 10^logP; %bar
        out = P * 14.5038; %psia
    end

t_prop_in = 308.706; %K (96 degF)
p_high_prop = prop_sat_pressure(319.817); %psia (116 degF)

t_prop_out = t_out + 2.8; % 5 degF above outlet air temrperature
p_low_prop = prop_sat_pressure(t_prop_out - 10.11); %psia (10.11 degC is 18.2 degF)

disp(t_out + 2.8)
% H in, 308.7060 K, 234.9027 psia
H_in_prop = 12.991; % kj/mol


% H out, 287.2849 K, 79.9789 psia
H_out_prop = 24.192;

% H after isentropic compression
H_comp_prop = 28.657;

MW_prop = 44.1; % g/mol

COP_calc(H_in_prop,H_out_prop,MW_prop,H_comp_prop)

disp(t_prop_in)
disp(p_high_prop)
disp(' ')
disp(t_prop_out)
disp(p_low_prop)

end
