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
