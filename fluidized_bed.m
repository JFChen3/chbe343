function fluidized_bed

%Collected data for water column
flow_rate_water = 0:0.1:1.2; %L/min
%mm H2O
pdrop_water = [0, 65, 110, 178, 226, 227, 230, 231, 232, 232, 232, 233, 233];
%mm
bed_height_water = [287, 287, 287, 287, 288, 298, 308, 317, 328, 334, 344, 354, 363];

%Collected data for air column
flow_rate_air = [0:20, 25]; %L/min at 1 bar abs, 20 C
%cm H2O
pdrop_air_lower = [0, 2.6, 6.0, 8.8, 12.0, 15.1, 18.3, 21.5, 24.7, 28.6, 32.2, 36.6, 39.2, 39.5, 40.7,41.5,41.5,41.5,41.5,41,41,41];
pdrop_air_upper = [0, 2.6, 6.0, 8.8, 12.0, 15.1, 18.3, 21.5, 24.7, 28.6, 32.2, 36.6, 39.2, 39.5, 40.7,42,42,42,42,42.5,42.5,42.5];
%mm
bed_height_air_lower = [316,316,316,316,316,316,316,316,316,316,316,316,316,317,318,323,328,330,335,340,340,360];
bed_height_air_upper = [316,316,316,316,316,316,316,316,316,316,316,316,316,317,318,323,328,335,340,350,360,390];

%Water column plot
plot(flow_rate_water,pdrop_water,'bo-',flow_rate_water,bed_height_water,'g^-')
legend('Pressure Drop','Bed Height','Location','Northwest')
axis([0,1.2,0,400])
xlabel('Flow Rate(L/min)')
ylabel('Pressure Drop(mm water), Bed Height(mm)')

%Air column plot
figure
hold on
plot(flow_rate_air,pdrop_air_lower,'bo-',flow_rate_air,bed_height_air_lower,'g^-')
legend('Pressure Drop','Bed Height','Location','Northwest')
plot(flow_rate_air,pdrop_air_upper,'bo-',flow_rate_air,bed_height_air_upper,'g^-')
xlabel('Flow Rate(L/min at 1 bar abs, 20^o C)')
ylabel('Pressure Drop(cm water), Bed Height(mm)')

%Calculate ballotini size
%Water
disp('Ballotini size, water column: ')
size_water = calc_ballotini_size(flow_rate_water, pdrop_water, bed_height_water, 1);
disp(size_water)

disp('Ballotini size, air')
size_air = calc_ballotini_size(flow_rate_air, pdrop_air_lower, bed_height_air_lower,2);
disp(size_air);
end

function [ballotini_size] = calc_ballotini_size(flow_rate, pdrop, bed_height, fluid)
% Inputs:   flow_rate in L/min
%           pdrop in mm H2O (water, fluid=1) or cm H2O (air, fluid=2)
%           bed_height in mm

column_diameter = 2.0*2.54/100; %column diameter in m
column_area = pi*(column_diameter/2)^2; %column cross sectional area, m^2

%Convert inputs to correct units
flow_rate = flow_rate/60000; %Convert L/min to m^3/s
bed_height = bed_height/1000; % Convert mm to m

bed_volume = bed_height.*column_area; %m^3

if fluid == 1
    pdrop = pdrop*9.80665; %Convert mm H2O to Pa
    epsilon = get_porosity(0.838,2.5e3,bed_volume);
    velocity = get_superficial_velocity(flow_rate, column_area); %m/s
    density = 1000; %kg/m^3
    viscosity = 1.002e-3;%Pa*s
else
    pdrop = pdrop*98.0665; %Convert cm H2O to Pa
    %======================== Ethan Added Part =====================================
    %Correct air flow rate to experimental conditions
    R = 8.314; % (L*kPa)/(K*mol) m3?Pa?K?1?mol?1
    P_original = 100; %kPa (1 bar)
    T_original = 293; %K
    T_experimental = 294.261; %K (70 degF)
    P_experimental = pdrop; %in cmH20
    P_experimental = 100 + P_experimental .* .0980665; % convert cmH20 guage to kPa absolute
    molar_flow_air = P_original .* flow_rate ./ (R * T_original);
    flow_rate_air_corrected = molar_flow_air .*  ( R * T_experimental ./ P_experimental );
    %====use flow_rate_air_corrected rather than flow_rate_air to use corrected values====
    velocity = get_superficial_velocity(flow_rate_air_corrected, column_area); %m/s
    P_atm = 101.325;
    density = P_atm*28.97/(R*T_experimental); %kg/m^3
    viscosity = 1.846e-5; %Pa*s
    epsilon = get_porosity(0.865,2.5e3,bed_volume);
end

ballotini_size = rearranged_ergun(bed_height,velocity,epsilon,density,viscosity,pdrop);
end

function [particle_diameter] = rearranged_ergun(L, Vsm, epsilon, rho, mu, h)
% Calculate ballotini size using rearranged Ergun equation for water column
% L = bed height (m)
% Vsm = superficial velocity (m/s)
% epsilon = porosity
% rho = density of fluid (kg/m^3)
% mu = viscosity of fluid (N*s/m^2)
% h = pressure drop (N/m^2)

a = -h;
b = (1.75*L.*Vsm.^2.*(1-epsilon).*rho)./epsilon.^3;
c = (150*L.*Vsm.*(1-epsilon).^2.*mu)./epsilon.^3;


particle_diameter = zeros(1,numel(a));
for i = 1:numel(a)-1
    solutions = roots([a(i+1) b(i+1) c(i+1)]);
    particle_diameter(i,:) = solutions(solutions >= 0); %Get the positive root
end

end

function [velocity] = get_superficial_velocity(flow_rate, column_area)

velocity = flow_rate/column_area;

end

function [epsilon] = get_porosity(mass_solid,density_solid,bed_volume)

epsilon = 1 - mass_solid./(density_solid.*bed_volume);

end