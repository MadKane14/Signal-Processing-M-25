function main 
clc;
clear;
close all;
syms t;

% Number of Fourier coefficients
N = 10; 

%% 1) x1(t) = 1 + sin(2*pi*f0*t) + 2*cos(2*pi*f0*t) + cos(4*pi*f0*t + pi/3)
f0 = 1; 
T1 = 1/f0; 
xt1 = 1 + sin(2*pi*f0*t) + 2*cos(2*pi*f0*t) + cos(4*pi*f0*t + pi/3);
X1 = fourierCoeff(t, xt1, T1, 0, T1, N);

%% 2) x2(t) = rect(t/T)  (rect defined over [-tau/2, tau/2])
T2 = 2; tau2 = 1; 
xt2 = piecewise(abs(t) <= tau2/2, 1, 0);
X2 = fourierCoeff(t, xt2, T2, -tau2/2, tau2/2, N);

%% 3) Triangular pulse (convolution of rect with itself)
T3 = 2; tau3 = 1;
xt3 = piecewise(-tau3 <= t & t < 0, 1 + t/tau3, ...
                 0 <= t & t < tau3, 1 - t/tau3, 0);
X3 = fourierCoeff(t, xt3, T3, -tau3, tau3, N);

%% 4) Impulse train sum_{n=-∞}^{∞} δ(t - nT)
T4 = 2; 
X4 = ones(1, 2*N + 1) / T4;  % Fourier coefficients for impulse train

% Represent xt4 in time domain for plotting: impulses at multiples of T
n_range = -3:3;                   
t_impulses = n_range * T4;           
amp_impulses = ones(size(t_impulses)); 

%% Plot Time-domain signals
t_vals = linspace(-2, 2, 400);
xt1_vals = double(subs(xt1, t, t_vals));
xt2_vals = double(subs(xt2, t, t_vals));
xt3_vals = double(subs(xt3, t, t_vals));

figure;
subplot(4,1,1), plot(t_vals, xt1_vals, 'LineWidth', 1.5), title('x₁(t)'), grid on
subplot(4,1,2), plot(t_vals, xt2_vals, 'LineWidth', 1.5), title('x₂(t)'), grid on
subplot(4,1,3), plot(t_vals, xt3_vals, 'LineWidth', 1.5), title('x₃(t)'), grid on
subplot(4,1,4), stem(t_impulses, amp_impulses, 'filled'), title('x₄(t) (Impulse Train)'), grid on

%% Plot Fourier-domain magnitude
k_vals = -N:N;
figure;
subplot(4,1,1), stem(k_vals, abs(X1), 'filled'), title('|X₁[k]|'), grid on
subplot(4,1,2), stem(k_vals, abs(X2), 'filled'), title('|X₂[k]|'), grid on
subplot(4,1,3), stem(k_vals, abs(X3), 'filled'), title('|X₃[k]|'), grid on
subplot(4,1,4), stem(k_vals, abs(X4), 'filled'), title('|X₄[k]|'), grid on

end

function X = fourierCoeff(t, xt, T, t1, t2, N)
% Compute Fourier series coefficients
w0 = 2*pi/T;
X = zeros(1, 2*N+1);

for k = -N:N
    xk = (1/T) * int(xt * exp(-1j*k*w0*t), t, t1, t2);
    X(k+N+1) = double(xk); 
end
end
