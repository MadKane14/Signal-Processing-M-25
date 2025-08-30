% Define the time vector 'n' for which we will compute x[n]
n = -10:10;

% --- Define Frequency-Domain Expressions ---
syms w; % Symbolic variable for function handles

% Part (a)
X_a_func = @(w) 1 + 3.*exp(-1j.*w) + 2.*exp(-2j.*w) - 4.*exp(-3j.*w) + exp(-10j.*w);

% Part (b)
X_b_func = @(w) cos(w).^2 + sin(3.*w).^2;

% Part (c) - NUMERICAL APPROXIMATION
% We cannot represent a true Dirac delta. We will approximate it with a
% very narrow and tall Gaussian function. This is expected to give an
% incorrect result, which is the point of the exercise for the report.
% A Gaussian with standard deviation 'sigma' has an area of sigma*sqrt(2*pi).
% To make the area 1 (like a delta), we must scale it by 1/(sigma*sqrt(2*pi)).
sigma = 0.001; % Make the Gaussian very narrow
gauss = @(w, mu) (1/(sigma*sqrt(2*pi))) * exp(-(w-mu).^2 / (2*sigma^2));

% Now, build the train of approximated deltas from -pi to pi
% The terms are for k = -2, -1, 0, 1, 2
X_c_func_approx = @(w) (-1)^(-2) * gauss(w, -pi) ...
                      + (-1)^(-1) * gauss(w, -pi/2) ...
                      + (-1)^(0)  * gauss(w, 0) ...
                      + (-1)^(1)  * gauss(w, pi/2) ...
                      + (-1)^(2)  * gauss(w, pi);

% Part (d)
X_d_func = @(w) (1 - (1/3).*exp(-1j.*w)) ./ (1 - (1/4).*exp(-1j.*w) - (1/8).*exp(-2j.*w));


% --- 1. Numerical Computation using the Inverse DTFT Integral ---
x_n_a_numerical = zeros(size(n));
x_n_b_numerical = zeros(size(n));
x_n_c_numerical = zeros(size(n)); % To store the incorrect result for (c)
x_n_d_numerical = zeros(size(n));

for i = 1:length(n)
    k = n(i);
    integrand_a = @(w) X_a_func(w) .* exp(1j.*w.*k);
    integrand_b = @(w) X_b_func(w) .* exp(1j.*w.*k);
    integrand_c = @(w) X_c_func_approx(w) .* exp(1j.*w.*k);
    integrand_d = @(w) X_d_func(w) .* exp(1j.*w.*k);
    
    x_n_a_numerical(i) = (1/(2*pi)) * integral(integrand_a, -pi, pi);
    x_n_b_numerical(i) = (1/(2*pi)) * integral(integrand_b, -pi, pi);
    x_n_c_numerical(i) = (1/(2*pi)) * integral(integrand_c, -pi, pi);
    x_n_d_numerical(i) = (1/(2*pi)) * integral(integrand_d, -pi, pi);
end

% Take the real part to remove small numerical errors
x_n_a_numerical = real(x_n_a_numerical);
x_n_b_numerical = real(x_n_b_numerical);
x_n_c_numerical = real(x_n_c_numerical);
x_n_d_numerical = real(x_n_d_numerical);


% --- 2. Analytical Computation ---
x_n_a_analytical = zeros(size(n));
x_n_a_analytical(n == 0) = 1; x_n_a_analytical(n == 1) = 3;
x_n_a_analytical(n == 2) = 2; x_n_a_analytical(n == 3) = -4;
x_n_a_analytical(n == 10) = 1;

x_n_b_analytical = zeros(size(n));
x_n_b_analytical(n == 0) = 1; x_n_b_analytical(n == 2) = 0.25;
x_n_b_analytical(n == -2) = 0.25; x_n_b_analytical(n == 6) = -0.25;
x_n_b_analytical(n == -6) = -0.25;

x_n_c_analytical = (1/(2*pi)) * (1 + cos(pi*n) - 2*cos(pi/2*n));
x_n_c_analytical = round(x_n_c_analytical, 4);

u_n = (n >= 0);
x_n_d_analytical = ((2/9)*(1/2).^n + (7/9)*(-1/4).^n) .* u_n;


% --- 3. Plotting the Results ---
% Plot for Numerically Computed Signals
figure('Name', 'Numerically Computed Inverse DTFTs', 'NumberTitle', 'off');
sgtitle('Time-Domain Signals (Numerical Integration)', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); stem(n, x_n_a_numerical, 'b', 'filled'); title('(a) Numerical'); grid on; axis tight;
subplot(2, 2, 2); stem(n, x_n_b_numerical, 'r', 'filled'); title('(b) Numerical'); grid on; axis tight;
subplot(2, 2, 3); stem(n, x_n_c_numerical, 'g', 'filled'); title('(c) Numerical (Incorrect)'); grid on; axis tight;
subplot(2, 2, 4); stem(n, x_n_d_numerical, 'm', 'filled'); title('(d) Numerical'); grid on; axis tight;

% Plot for Analytically Derived Signals
figure('Name', 'Analytically Derived Inverse DTFTs', 'NumberTitle', 'off');
sgtitle('Time-Domain Signals (Analytical Formulas)', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); stem(n, x_n_a_analytical, 'b', 'filled'); title('(a) Analytical'); grid on; axis tight;
subplot(2, 2, 2); stem(n, x_n_b_analytical, 'r', 'filled'); title('(b) Analytical'); grid on; axis tight;
subplot(2, 2, 3); stem(n, x_n_c_analytical, 'g', 'filled'); title('(c) Analytical (Correct)'); grid on; axis tight;
subplot(2, 2, 4); stem(n, x_n_d_analytical, 'm', 'filled'); title('(d) Analytical'); grid on; axis tight;
