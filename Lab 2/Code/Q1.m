f = linspace(-0.5, 0.5, 1000);

% --- Signal (a): Rectangular Pulse ---
n_a = -2:10;
x_a = (n_a >= 2) - (n_a >= 6); % x[n] is 1 for n=2,3,4,5
X_a = dtft(n_a, x_a, f);

% --- Signal (b): Left-Sided Exponential with |n| 
n_b = -20:-2; 
x_b = (1/3).^(abs(n_b));
X_b = dtft(n_b, x_b, f);

% --- Signal (c): Finite Ramp ---
n_c = -3:3;
x_c = n_c;
X_c = dtft(n_c, x_c, f);

% --- Signal (d): Modulated Sinc Function ---
n_d = -50:50; 
x_d = sin(pi*n_d/5) ./ (pi*n_d);
% Handle the n=0 case (limit of sin(ax)/(ax) is 1, so sin(ax)/(pi*x) -> a/pi)
x_d(n_d == 0) = (1/5); 
% Apply cosine modulation
x_d = x_d .* cos(pi*n_d/2);
X_d = dtft(n_d, x_d, f);


% --- Plotting Section ---

% Figure 1: Time-Domain Signals
figure('Name', 'Time-Domain Signals');
subplot(2,2,1);
stem(n_a, x_a, 'filled');
title('(a) x[n] = u[n-2] - u[n-6]');
xlabel('n (sample)'); ylabel('x[n]');
grid on; axis([-2 10 -0.5 1.5]);

subplot(2,2,2);
stem(n_b, x_b, 'filled');
title('(b) x[n] = (1/3)^{|n|} u[-n-2]'); 
xlabel('n (sample)'); ylabel('x[n]');
grid on;

subplot(2,2,3);
stem(n_c, x_c, 'filled');
title('(c) x[n] = n, for -3 \leq n \leq 3');
xlabel('n (sample)'); ylabel('x[n]');
grid on;

subplot(2,2,4);
stem(n_d, x_d, 'filled');
title('(d) x[n] = sinc(n/5)cos(\pi n/2)');
xlabel('n (sample)'); ylabel('x[n]');
grid on; axis([-20 20 -0.2 0.3]);

% Figure 2: Magnitude Spectra
figure('Name', 'Magnitude Spectra');
subplot(2,2,1);
plot(f, abs(X_a));
title('(a) Magnitude |X(e^{j\omega})|');
xlabel('Normalized Frequency (f)'); ylabel('Magnitude');
grid on;

subplot(2,2,2);
plot(f, abs(X_b));
title('(b) Magnitude |X(e^{j\omega})|'); 
xlabel('Normalized Frequency (f)'); ylabel('Magnitude');
grid on;

subplot(2,2,3);
plot(f, abs(X_c));
title('(c) Magnitude |X(e^{j\omega})|');
xlabel('Normalized Frequency (f)'); ylabel('Magnitude');
grid on;

subplot(2,2,4);
plot(f, abs(X_d));
title('(d) Magnitude |X(e^{j\omega})|');
xlabel('Normalized Frequency (f)'); ylabel('Magnitude');
grid on; axis([-0.5 0.5 -0.1 0.6]);

% Figure 3: Phase Spectra
figure('Name', 'Phase Spectra');
subplot(2,2,1);
plot(f, angle(X_a));
title('(a) Phase \angleX(e^{j\omega})');
xlabel('Normalized Frequency (f)'); ylabel('Phase (radians)');
grid on;

subplot(2,2,2);
plot(f, angle(X_b));
title('(b) Phase \angleX(e^{j\omega})'); 
xlabel('Normalized Frequency (f)'); ylabel('Phase (radians)');
grid on;

subplot(2,2,3);
plot(f, angle(X_c));
title('(c) Phase \angleX(e^{j\omega})');
xlabel('Normalized Frequency (f)'); ylabel('Phase (radians)');
grid on;

subplot(2,2,4);
plot(f, angle(X_d));
title('(d) Phase \angleX(e^{j\omega})');
xlabel('Normalized Frequency (f)'); ylabel('Phase (radians)');
grid on; axis([-0.5 0.5 -pi pi]);


% --- DTFT Function Definition ---
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