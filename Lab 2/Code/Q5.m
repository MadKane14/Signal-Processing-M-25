%% Part (b): Generate Sine Wave and its Noisy Version

% --- Signal Parameters ---
N = 200;         % Number of samples
w0 = pi/10;      % Angular frequency of the sine wave
n = 0:N-1;       % Time index vector

% --- Generate Signals ---
s = 5 * sin(w0 * n);          % Original sine wave
w = randn(1, N);              % Gaussian white noise (mean=0, variance=1)
x = s + w;                    % Noisy signal

% --- Plotting ---
figure('Name', 'Time-Domain Analysis', 'NumberTitle', 'off');

% First tile: Plot original and noisy signals
subplot(2, 2, 1);
plot(n, s, 'b-', 'LineWidth', 1.5);
hold on;
plot(n, x, 'r-', 'LineWidth', 0.8);
hold off;
title('Original vs. Noisy Signal');
xlabel('Time index (n)');
ylabel('Amplitude');
legend('s[n] (Original)', 'x[n] (Noisy)');
grid on;
axis tight;

%% Part (c): Filter the Signal with Moving Average Filter

M_values = [5, 21, 51]; % Filter orders
denoised_signals = cell(1, length(M_values)); % Store results

for i = 1:length(M_values)
    M = M_values(i);
    
    % Define the impulse response (the filter)
    h = (1/M) * ones(1, M);
    
    % Filter the noisy signal using convolution
    % Using 'full' as requested. The output length will be N + M - 1.
    y = conv(x, h, 'full');
    
    % Store the first N samples for plotting and DTFT analysis
    denoised_signals{i} = y(1:N);
    
    % Plot in the remaining 3 tiles
    subplot(2, 2, i + 1);
    plot(n, s, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(n, denoised_signals{i}, 'g-', 'LineWidth', 1.2);
    hold off;
    title(['Denoised Signal with M = ', num2str(M)]);
    xlabel('Time index (n)');
    ylabel('Amplitude');
    legend('s[n] (Original)', ['y[n] (Denoised, M=', num2str(M),')']);
    grid on;
    axis tight;
end

%% Part (e): DTFT Analysis

% --- Frequency Vector ---
% Define a high-resolution frequency vector for smooth plots
f = -0.5:0.001:0.5;

% --- Create a new figure for frequency-domain plots ---
figure('Name', 'Frequency-Domain (DTFT) Analysis', 'NumberTitle', 'off');

% --- DTFT of the Noisy Signal ---
X = dtft(n, x, f);
subplot(2, 2, 1);
plot(f, abs(X), 'r-');
title('DTFT of Noisy Signal x[n]');
xlabel('Normalized Frequency (f)');
ylabel('|X(f)|');
grid on;
axis tight;

% --- DTFT of the Denoised Signals ---
for i = 1:length(M_values)
    M = M_values(i);
    y_denoised = denoised_signals{i};
    
    % Calculate DTFT
    Y = dtft(n, y_denoised, f);
    
    % Plot
    subplot(2, 2, i + 1);
    plot(f, abs(Y), 'g-');
    title(['DTFT of Denoised Signal, M = ', num2str(M)]);
    xlabel('Normalized Frequency (f)');
    ylabel('|Y(f)|');
    grid on;
    axis tight;
end

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
