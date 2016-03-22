function flow_meters

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