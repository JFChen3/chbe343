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
