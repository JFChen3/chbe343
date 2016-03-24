function flow_meters
close all

%Collected data
volumes = [21, 21, 21, 11, 11, 5, 2, 2]; %L
times = [47, 32, 40, 37, 55, 43, 31, 63]; %s
%in H2O
diff_head_o = [7.1, 15.4, 10.2, 4.5, 1.98, 0.49, 0.10, 0.02];
head_loss_o = [2.7, 6.3, 4.0, 1.77, 0.78, 0.19, 0.01, 0.01];
diff_head_v = [21.0, 40.4, 27.3, 14.5, 6.0, 1.92, 0.32, 0.1];
head_loss_v = [3.3, 5.7, 4.25, 2.8, 1.8, 0.73, 0.15, 0.05];
diff_head_p = [3.3, 6.5, 5.1, 2.4, 1.2, 0.3, 0.1, 0.03];
%gal/min
coriolis = [8.1, 11.3, 9.2, 6.3, 4.3, 2.3, 1.0, 0.5];

%%Calculate flow rates, m^3/s
flow_rate_o = orifice(0.62, 20e-3, 24e-3, diff_head_o*2.54e-2); %Orifice
flow_rate_v = orifice(0.98, 14e-3, 24e-3, diff_head_v*2.54e-2); %Venturi
flow_rate_p = pitot(24e-3, diff_head_p*2.54e-2); %Pitot
coriolis_nice_units = coriolis./15850.3;
measured = volumes./(1000*times);
disp('Flow rates: Measured, Coriolis, Orifice Plate, Venturi, Pitot Tube')
disp([measured' coriolis_nice_units' flow_rate_o' flow_rate_v' flow_rate_p'])

%%Calculate error from measured flow rate
error_orifice = abs(flow_rate_o - measured)./measured;
error_venturi = abs(flow_rate_v - measured)./measured;
error_pitot = abs(flow_rate_p - measured)./measured;
error_coriolis = abs(coriolis_nice_units - measured)./measured;
disp('Relative Error to Measured Flow Rate: Orifice Plate, Venturi, Pitot Tube')
disp([error_orifice', error_venturi', error_pitot'])
disp('Relative Error to Measured Flow Rate: Coriolis Meter')
disp(error_coriolis')


% Part 4: compare head loss for orifice and venturi
measured_plot = reorder(measured);
head_loss_o_plot = reorder(head_loss_o);
head_loss_v_plot = reorder(head_loss_v);
diff_head_o_plot = reorder(diff_head_o);
diff_head_v_plot = reorder(diff_head_v);

figure
hold on
plot(measured_plot, head_loss_o_plot, 'r--', 'LineWidth', 1.2);
plot(measured_plot, head_loss_v_plot, 'g-x', 'LineWidth', 1.2)
xlabel('Flow Rate (m^3/s)')
ylabel('Head Loss (inches H2O)')
legend('Orifice Plate', 'Venturi Meter', 'Location', 'SouthEast')

avg_head_loss_v = sum(head_loss_v)/length(head_loss_v);
avg_head_loss_o = sum(head_loss_o)/length(head_loss_o);

fprintf('Average head loss for Venturi meter:\t %.2f in H2O\n', avg_head_loss_v)
fprintf('Average head loss for orifice meter:\t %.2f in H2O\n', avg_head_loss_o)

%Part 5: compare differential head for venturi and oriface plate
figure
hold on
plot(measured_plot, diff_head_o_plot, 'r--', 'LineWidth', 1.2)
plot(measured_plot, diff_head_v_plot, 'g-x', 'LineWidth', 1.2)
xlabel('Flow Rate (m^3/s)')
ylabel('Differential Head (inches H2O)')
legend('Orifice Plate', 'Venturi Meter', 'Location', 'NorthWest')

avgscale = mean(diff_head_v_plot./diff_head_o_plot);

end

function [flow_rate] = orifice(Cd, d0, d1, head)

A0 = pi*(d0/2)^2;
A1 = pi*(d1/2)^2;
flow_rate = Cd.*A0.*((1-(A0./A1).^2).^(1/2)).*sqrt(2.*9.81.*head);

end

function [flow_rate] = pitot(d, head)

A = pi*(d/2)^2;
flow_rate = A.*sqrt(2.*9.81.*head);

end

function [ordered] = reorder(vector)
%reorders the vectors for plotting

ordered = [vector(2:3) vector(1) vector(4:end)];
end

rho = 998.02;
mu = 1.002e-3;
dia = .24; %m


flow_rates = [
    0.4468    0.5110    0.2636    0.4589    0.5801
    0.6563    0.7129    0.3883    0.6365    0.8142
    0.5250    0.5804    0.3160    0.5232    0.7212
    0.2973    0.3975    0.2099    0.3813    0.4947
    0.2000    0.2713    0.1392    0.2453    0.3498
    0.1163    0.1451    0.0693    0.1388    0.1749
    0.0645    0.0631    0.0313    0.0566    0.1010
    0.0317    0.0315    0.0140    0.0317    0.0553   
    ];


flow_rates = sort(flow_rates,1);
flow_rates = flow_rates .* 10^-3;

velocities = flow_rates ./ ((dia/2)^2 * pi);

re = rho.*velocities.*dia./mu;

v_turb = mu*4000/(rho*dia);
v_lam = mu*2000/(rho*dia);

q_turb = v_turb * (dia/2)^2 * pi;
q_lam = v_lam * (dia/2)^2 * pi;

coriolis = flow_rates(:,2);
stopwatch = flow_rates(:,1);
orifice = flow_rates(:,3);
venturi = flow_rates(:,4);
pitot = flow_rates(:,5);

coriolis_v = velocities(:,2);
stopwatch_v = velocities(:,1);
orifice_v = velocities(:,3);
venturi_v = velocities(:,4);
pitot_v = velocities(:,5);

coriolis_error = (coriolis - stopwatch) ./ stopwatch;
orifice_error = (orifice - stopwatch) ./ stopwatch;
venturi_error = (venturi - stopwatch) ./ stopwatch;
pitot_error = (pitot - stopwatch) ./ stopwatch;

stopwatch_error_c = (stopwatch - coriolis) ./ coriolis;
orifice_error_c = (orifice - coriolis) ./ coriolis;
venturi_error_c = (venturi - coriolis) ./ coriolis;
pitot_error_c = (pitot - coriolis) ./ coriolis;

mean_coriolis_error = mean(coriolis_error);
mean_orifice_error = mean(orifice_error);
mean_venturi_error = mean(venturi_error);
mean_pitot_error = mean(pitot_error);

%disp(mean_coriolis_error)
%disp(mean_orifice_error)
%disp(mean_venturi_error)
%disp(mean_pitot_error)

% pitot plot
y_range = -80:80;

hold on
%plot(stopwatch,coriolis_error.*100)
plot(stopwatch,orifice_error.*100,'--r')
plot(stopwatch,venturi_error.*100,'-xg')
%plot(stopwatch,pitot_error.*100,'-sk')
plot(ones(1,numel(y_range)).*q_turb,y_range,'-.m')
plot(ones(1,numel(y_range)).*q_lam,y_range,'-.m')
%title('Error as a Function of Flow Rate')
xlabel('Flow Rate (m^3/s)')
ylabel('Percent Error Compared to Stopwatch Measurement')
%legend('Coriolis Meter','Orifice Plate','Venturi Meter','Pitot Tube')
legend('Orifice Plate','Venturi Meter','Bounds of Transition Regime','location','northwest')

figure()
% pitot plot
hold on
%plot(coriolis,stopwatch_error_c.*100)
plot(coriolis,orifice_error_c.*100,'--r')
plot(coriolis,venturi_error_c.*100,'-xg')
%plot(coriolis,pitot_error_c.*100,'-sk')
plot(ones(1,numel(y_range)).*q_turb,y_range,'-.m')
plot(ones(1,numel(y_range)).*q_lam,y_range,'-.m')
%title('Error as a Function of Flow Rate')
xlabel('Flow Rate (m^3/s)')
ylabel('Percent Error Compared to Coriolis Measurement')
%legend('Stopwatch Measurement','Orifice Plate','Venturi Meter','Pitot Tube')
legend('Orifice Plate','Venturi Meter','Bounds of Transition Regime','location','northwest')
