function main
syms t;

% Number of harmonics
N = 10;

f0 = 1; 
T1 = 1/f0;
xt1 = 1 + sin(2*pi*f0*t) + 2*cos(2*pi*f0*t) + cos(4*pi*f0*t + pi/4);
X1 = fourierCoeff(t, xt1, T1, 0, T1, N);

T2 = 2; tau = 1;
xt2 = piecewise(abs(t) <= tau/2, 1, 0);
X2 = fourierCoeff(t, xt2, T2, -tau/2, tau/2, N);

T3 = 2; tau = 1;
xt3 = piecewise(-tau <= t & t < 0, 1 + t/tau, ...
                 0 <= t & t < tau, 1 - t/tau, 0);
X3 = fourierCoeff(t, xt3, T3, -tau, tau, N);

T4 = 2;  
X4 = ones(1, 2*N + 1)/T4;

time1 = linspace(0, T1, 500);
time2 = linspace(-T2/2, T2/2, 500);
time3 = linspace(-T3/2, T3/2, 500);
time4 = linspace(-T4, T4, 500); 

xhat1 = partialFourierSum(X1, T1, time1);
xhat2 = partialFourierSum(X2, T2, time2);
xhat3 = partialFourierSum(X3, T3, time3);
xhat4 = partialFourierSum(X4, T4, time4);

xorig1 = double(subs(xt1, t, time1));
xorig2 = double(subs(xt2, t, time2));
xorig3 = double(subs(xt3, t, time3));

xorig4 = zeros(size(time4));
impulse_positions = -T4:T4:T4; 
for pos = impulse_positions
    [~, idx] = min(abs(time4 - pos));
    xorig4(idx) = 1; 
end

figure;

subplot(2,2,1);
plot(time1, xorig1, 'LineWidth', 1.5); hold on;
plot(time1, xhat1, '--', 'LineWidth', 1.5);
title('Signal 1'); legend('Original', 'Reconstructed');

subplot(2,2,2);
plot(time2, xorig2, 'LineWidth', 1.5); hold on;
plot(time2, xhat2, '--', 'LineWidth', 1.5);
title('Signal 2'); legend('Original', 'Reconstructed');

subplot(2,2,3);
plot(time3, xorig3, 'LineWidth', 1.5); hold on;
plot(time3, xhat3, '--', 'LineWidth', 1.5);
title('Signal 3'); legend('Original', 'Reconstructed');

subplot(2,2,4);
stem(time4, xorig4, 'LineWidth', 1.5); hold on;
plot(time4, xhat4, '--', 'LineWidth', 1.5);
title('Impulse Train'); legend('Original', 'Reconstructed');

end

function X = fourierCoeff(t, xt, T, t1, t2, N)
% function to compute Fourier series coefficients
% t   - symbolic variable for time
% xt  - continuous-time signal (symbolic expression)
% T   - period of the signal
% t1,t2 - valid integration limits
% N   - coefficients will be computed for k = -N : N
% X   - (2N+1) length Fourier series coefficients vector

w0 = 2*pi/T;

X = zeros(1, 2*N+1);

for k = -N:N
    xk = (1/T) * int(xt * exp(-1j*k*w0*t), t, t1, t2);
    X(k+N+1) = double(xk); 
end
end

function x_hat = partialFourierSum(A, T, time_grid)
% partialFourierSum - Performs partial Fourier reconstruction
%
% Inputs:
%   A         - (2N+1)-length vector of Fourier coefficients C_k
%   T         - Period of signal
%   time_grid - Vector of time samples for reconstruction
%
% Output:
%   x_hat     - Reconstructed signal values at time_grid

N = (length(A)-1)/2;

w0 = 2*pi/T;

x_hat = zeros(size(time_grid));

for k = -N:N
    Ck = A(k+N+1);        
    x_hat = x_hat + Ck * exp(1j*k*w0*time_grid);
end

x_hat = real(x_hat);
end