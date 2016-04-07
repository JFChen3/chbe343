function kinetics

%data for setting A, T=29.5 C
time_A = 20:20:360;
abs_A = [0.462, 0.337, 0.239, 0.169, 0.120, 0.085, 0.061, 0.043, 0.031, 0.023, 0.017, 0.012, 0.009, 0.007, 0.005, 0.004, 0.003, 0.003];
[k1star_A, k2_A] = compute_all(time_A, abs_A, 0.5, 29.5);

%data for setting B, T=36.4 C
time_B = 20:20:320;
abs_B = [0.520, 0.308, 0.181, 0.110, 0.062, 0.036, 0.021, 0.011, 0.005, 0.001, 0.000, -0.004, -0.006, -0.007, -0.008, -0.008];
[k1star_B, k2_B] = compute_all(time_B(1:numel(time_B)-6), abs_B(1:numel(time_B)-6), 8, 36.4);

%data for setting C, T=48.5 C
time_C = 15:5:90;
abs_C = [0.001, 0.053, 0.111, 0.085, 0.060, 0.043, 0.032, 0.025, 0.021, 0.020, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018];
[k1star_C, k2_C] = compute_all(time_C(3:numel(time_C)-5), abs_C(3:numel(time_C)-5), 0.2, 48.5);

%%Activation Energies
[E1, uncert_1, preexp1, exp_uncert1] = calc_Ea([k1star_A, k1star_B, k1star_C],[29.5, 36.4, 48.5]+273.15);
[E2, uncert_2, preexp2, exp_uncert2] = calc_Ea([k2_A, k2_B, k2_C],[29.5, 36.4, 48.5]+273.15);
Hrxn = E1-E2;

fprintf('-----------------------------------\n')
fprintf('Activation Energies (J/mol) \n')
fprintf('-----------------------------------\n')
fprintf('E1                       %4.4e \n', E1)
fprintf('Uncertainty              %4.4e \n', uncert_1)
fprintf('Pre-Exponential          %4.4e \n', preexp1)
fprintf('Uncertainty              %4.4e \n', exp_uncert1)
fprintf('-----------------------------------\n')
fprintf('E2                       %4.4e \n', E2)
fprintf('Uncertainty              %4.4e \n', uncert_2)
fprintf('Pre-Exponential          %4.4e \n', preexp2)
fprintf('Uncertainty              %4.4e \n', exp_uncert2)
fprintf('-----------------------------------\n')
fprintf('Heat of Reaction: %4.4e J/mol \n', Hrxn)

end

function concentration = conc_from_rate(k1,k2,c10,t)
    concentration = k2/(k1+k2)*c10 + k1/(k1+k2).*exp(-(k1+k2).*t).*c10;
end

function [k1star, k2] = compute_all(time, abs, c_factor, T)

[A0, Aeq, slope] = calc_A0(time, abs, c_factor, T);
[k1star, k2] = rxn_constants(slope, A0, Aeq);

%Uncertainty
[sigma_m, sigma_b] = linear_reg_uncertainty(time,log((abs-Aeq)/(A0-Aeq)),slope,0);

%Check c2 at equilibrium
c2_eq = check_c2(A0, Aeq, k1star, k2);

fprintf('------------------------------\n')
fprintf('Results for T=%.1f C \n',T)
fprintf('------------------------------\n')
fprintf('A0                     %.4f\n', A0)
fprintf('Aeq                    %.4f\n', Aeq)
fprintf('k1*                    %4.4e\n', k1star)
fprintf('k2                     %4.4e\n', k2)
fprintf('Uncertainty, slope     %4.4e\n', sigma_m)
fprintf('Uncertainty, intercept %4.4e\n', sigma_b)
fprintf('Equilibrium c2         %.4f\n', c2_eq)

figure
plot(time, abs, '.-')
xlabel('Time (s)')
ylabel('Absorbance')
legend(sprintf('T=%.1f C',T))

% Plot to validate the kinetic model
conc = conc_from_rate(k1star,k2,A0,time);

figure
hold on
plot(time,abs,'rx')
plot(time,conc)
xlabel('Time (s)')
ylabel('Absorbance')
legend(sprintf('T=%.1f C, Data',T),'Model Generated from Rate Constants')

end

function [A0, Aeq, slope] = calc_A0(time, A, c, T)

Aeq = A(end) - c*A(end); %Estimate equilibrium absorbance value

tol = 0.0001;
error = 1;
count = 0;

A0_guess = max(A); %Guess initial A0

while error > tol
    log_A = log((A-Aeq)/(A0_guess-Aeq));
    p = polyfit(time, log_A, 1);
    A0_guess = A0_guess + p(2);
    error = abs(p(2));
    count = count+1;
    if count > 100
        disp('Took too many iterations to find A0, check inputs')
        break
    end
end

A0 = A0_guess;
slope = p(1);

figure
hold on
plot(time, log_A, 'o')
plot(time, p(1).*time + p(2))
xlabel('Time (s)')
ylabel('log((A-Aeq)/(A0-Aeq))')
legend(sprintf('Data, T=%.1f C',T), sprintf('Fit, y=%.4f*t',slope))

end

function [k1_star, k2] = rxn_constants(slope, A0, Aeq)
%Calculate rate constants
C = (A0 - Aeq)/Aeq;
k2 = -slope/(1+C);
k1_star = -slope-k2;
end

function [c2_eq] = check_c2(A0, Aeq, k1, k2)

K = k1/k2;
c2_eq = (A0-Aeq)/(K*Aeq);

end

function [Ea, uncert, preexp, exp_uncert] = calc_Ea(k, T)
%Temp must be in Kelvin

x = 1./(8.314*T);
y = log(abs(k));

p = polyfit(x, y, 1);
Ea = -p(1);
preexp = exp(p(2));

figure
plot(x,y,'o',x,p(1)*x+p(2),'-')
xlabel('1/RT')
ylabel('ln(|k|)')
legend('Data',sprintf('Fit, y=%.1fx+%.1f',p(1),p(2)))

[uncert, exp_uncert] = linear_reg_uncertainty(x, y, p(1), p(2));

exp_uncert = abs(preexp*exp_uncert);

end

function [sigma_m, sigma_b] = linear_reg_uncertainty(x_arr,y_arr,m,b)
% Outputs the uncertainty in the m and b constants of y = mx + b
% x_arr, and y_arr are the the arrays of x and y experimental values
% m and b are constants from the linear interpolation

N = numel(x_arr);

x_sq = x_arr.^2;

delta = N .* sum(x_sq) - sum(x_arr).^2;

% Uncertainty in the measurements of y
sigma_y = sqrt(1/(N-2) .* sum((y_arr - b - m.*x_arr).^2));

% Uncertainty in b (y = mx + b)
sigma_b = sigma_y .* sqrt(sum(x_sq)./delta);

% Uncertainty in m (y = mx+b)
sigma_m = sigma_y .* sqrt(N./delta);

end
