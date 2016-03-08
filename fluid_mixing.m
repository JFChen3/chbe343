function fluid_mixing

%Given lever arm distance (m)
lever_arm = 0.11;

% Rotational speed in rpm, force in N
% Impeller size 62x65 mm
rpm_6265 = [64.5, 114.5, 134.0, 152.1, 170.5, 195.0, 212.0]; 
force_6265 = [0.45, 2.0, 3.8, 5.3, 6.0, 7.9, 9.9];
[rot_speed_6265, torque_6265, power_6265] = calc_params(rpm_6265, force_6265, lever_arm);
disp('---------------------------------------------------------------------------')
disp('Calculated Parameters, 62x65mm impeller')
disp('---------------------------------------------------------------------------')
disp('Rotation, rpm   Force, N   Rotational Speed, rad/s   Torque, N*m   Power, W')
for i = 1:numel(rpm_6265)
    fprintf('%.1f             %.1f       %.1f                      %.1f          %.1f \n',rpm_6265(i), force_6265(i), rot_speed_6265(i), torque_6265(i), power_6265(i))
end
disp('---------------------------------------------------------------------------')

% Impeller size 62x43 mm
rpm_6243 = [46.3, 95.4, 130.7, 153.9, 174.9, 199.5, 235.4, 250.0];
force_6243 = [1.1, 1.9, 3.1, 4.2, 4.8, 6.3, 8.6, 9.6];
[rot_speed_6243, torque_6243, power_6243] = calc_params(rpm_6243, force_6243, lever_arm);
disp('---------------------------------------------------------------------------')
disp('Calculated Parameters, 62x43mm impeller')
disp('---------------------------------------------------------------------------')
disp('Rotation, rpm   Force, N   Rotational Speed, rad/s   Torque, N*m   Power, W')
for i = 1:numel(rpm_6243)
    fprintf('%.1f             %.1f       %.1f                      %.1f          %.1f \n',rpm_6243(i), force_6243(i), rot_speed_6243(i), torque_6243(i), power_6243(i))
end
disp('---------------------------------------------------------------------------')

% Impeller size 62x21 mm
rpm_6221 = [72.0, 93.1, 128.0, 158.5, 183.0, 208.0, 252.0, 287.0, 309.4, 338.0, 367.1, 390.3];
force_6221 = [0.7, 0.7, 1.7, 2.4, 2.6, 3.4, 4.7, 5.3, 6.3, 7.4, 8.5, 9.4];
[rot_speed_6221, torque_6221, power_6221] = calc_params(rpm_6221, force_6221, lever_arm);
disp('---------------------------------------------------------------------------')
disp('Calculated Parameters, 62x21mm impeller')
disp('---------------------------------------------------------------------------')
disp('Rotation, rpm   Force, N   Rotational Speed, rad/s   Torque, N*m   Power, W')
for i = 1:numel(rpm_6221)
    fprintf('%.1f             %.1f       %.1f                      %.1f          %.1f \n',rpm_6221(i), force_6221(i), rot_speed_6221(i), torque_6221(i), power_6221(i))
end
disp('---------------------------------------------------------------------------')

%Power vs speed plot
hold on
plot(rot_speed_6265, power_6265, 'ro-')
plot(rot_speed_6243, power_6243, 'g^-')
plot(rot_speed_6221, power_6221, 'b.-')
xlabel('Rotational speed, rad/s')
ylabel('Power, W')
h = legend('62 x 65','62 x 43','62 x 21','Location','Southeast');
v = get(h,'title');
set(v,'string','Impeller Size, mm');

%Power Number vs Reynolds Number
density = 850; %mineral oil, kg/m^3
viscosity = 0.0153; %Pa*s, DOUBLE CHECK THIS VALUE
% all_rps = 60*[rpm_6265 rpm_6243 rpm_6221];
% all_power = [power_6265 power_6243 power_6221];
[Re_6265, P_6265] = dimensionless_groups(lever_arm*2, 60*rpm_6265, density, viscosity, power_6265);
[Re_6243, P_6243] = dimensionless_groups(lever_arm*2, 60*rpm_6243, density, viscosity, power_6243);
[Re_6221, P_6221] = dimensionless_groups(lever_arm*2, 60*rpm_6221, density, viscosity, power_6221);

% [reynolds_number, power_number] = dimensionless_groups(lever_arm*2, all_rps, density, viscosity, all_power);
% reynolds_number
% power_number
figure
plot(Re_6265, P_6265, 'r.-', Re_6243, P_6243, 'g.-', Re_6221, P_6221, 'b.-')
xlabel('Reynolds Number')
ylabel('Power Number')

end

function [rot_speed, torque, power] = calc_params(rpm, force, lever_arm)
% Calculate rotational speed, torque, and power

rot_speed = 60.*rpm./(2*pi); %rad/s
torque = force.*lever_arm; %N*m
power = 2*pi*torque.*rot_speed; %W

end

function [reynolds_number, power_number] = dimensionless_groups(D_impeller, N_rps, density, viscosity, power)

reynolds_number = ((D_impeller.^2).*N_rps.*density)./viscosity;
power_number = power./((N_rps.^3) .* (D_impeller.^5) .* density);

end