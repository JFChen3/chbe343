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

end

function [particle_diameter] = calc_ballotini_size()

end
function [velocity] = get_superficial_velocity(flow_rate, column_diameter)

column_area = pi*(column_diameter/2)^2;
velocity = flow_rate/column_area;

end

function [epsilon] = get_porosity(mass_solid,density_solid,bed_volume)

epsilon = 1 - mass_solid/(density_solid*bed_volume);

end