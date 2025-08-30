% Signal Used: x[n] = u[n-2] - u[n-6]
% This is a rectangular pulse of length 4, starting at n=2.

clear; clc; close all;

% --- 1. Define the Signal and Time/Frequency Vectors ---

% Define a discrete time vector 'n' over a sufficient range.
n = -10:10;

% Define the unit step function u[n]
u = @(n_vec) (n_vec >= 0);

% Define the signal x[n] = u[n-2] - u[n-6]
x = u(n-2) - u(n-6);

% Define the frequency vector for the DTFT calculation.
% We'll compute the DTFT at 1024 points between -0.5 and 0.5 (normalized frequency)
num_points = 1024;
f = linspace(-0.5, 0.5, num_points);

% --- 2. Method 1: PSD from Autocorrelation ---

% (a) Compute the autocorrelation of the signal using xcorr.
% The 'lags' vector from xcorr will serve as the time index for Rxx.
[Rxx, lags] = xcorr(x);

% (b) Find the PSD by computing the DTFT of the autocorrelation function
% using the custom 'dtft' function.
Sx_from_Rxx_complex = dtft(lags, Rxx, f);
Sx_from_Rxx = abs(Sx_from_Rxx_complex); % We only need the magnitude


% --- 3. Method 2: PSD from Signal's DTFT ---

% (a) Find the DTFT of the signal x[n] using the custom 'dtft' function.
X_omega_complex = dtft(n, x, f);

% (b) Find the PSD using the formula Sx(w) = |X(w)|^2.
Sx_from_X = abs(X_omega_complex).^2;


% --- 4. Plot the Results for Verification ---

% Define the angular frequency axis for plotting, from -pi to pi.
omega = 2 * pi * f;

% Create a 2x1 plot figure
figure('Name', 'PSD Verification');

% Plot PSD from Autocorrelation
subplot(2, 1, 1);
plot(omega, Sx_from_Rxx, 'g', 'LineWidth', 1.5);
title('Method 1: PSD from DTFT of Autocorrelation R_{xx}[n]');
xlabel('Frequency (\omega)');
ylabel('S_x(\omega)');
grid on;
xlim([-pi pi]);

% Plot PSD from Signal's DTFT
subplot(2, 1, 2);
plot(omega, Sx_from_X, 'r', 'LineWidth', 1.5);
title('Method 2: PSD from Squared Magnitude |X(\omega)|^2');
xlabel('Frequency (\omega)');
ylabel('S_x(\omega)');
grid on;
xlim([-pi pi]);

function X = dtft(n, x, f)
    % Evaluates the discrete-time Fourier transform of a given signal
    % n - vector of time indices where the signal is defined
    % x - discrete-time signal vector
    % f - vector of frequency values where DTFT is computed
    % X - complex vector of Fourier domain values
    
    % Ensure f is a column vector for matrix multiplication
    f = f(:);
    
    % Create a matrix of complex exponentials
    % Each row corresponds to a frequency in f
    % Each column corresponds to a time index in n
    exp_matrix = exp(-1j * 2 * pi * f * n);
    
    % Compute the DTFT via matrix-vector multiplication
    X = exp_matrix * x(:);
end