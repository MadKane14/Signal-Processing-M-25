function main
    % Parameters
    T = 2;              
    tau = 1;            
    t = linspace(-T/2, T/2, 2000);  

    % Original periodic rect function
    x_t = double(abs(t) <= tau/2);
    
    N_values = 1:100; 
    e_MAE = zeros(size(N_values));

    for idx = 1:length(N_values)
        N = N_values(idx);

        % Compute Fourier coefficients (numeric)
        Ck = fourierCoeff(t, x_t, T, N);

        % Reconstruct using partial sum
        x_hat = partialFourierSum(Ck, T, t);

        % Numerical integration for MAE
        e_MAE(idx) = (1/T) * trapz(t, abs(x_t - x_hat));
    end

    % Plot
    figure;
    plot(N_values, e_MAE, 'LineWidth', 1.5);
    xlabel('N (Number of Harmonics)');
    ylabel('Mean Absolute Error (MAE)');
    title('MAE of Rectangular Pulse Reconstruction vs N');
    grid on;
end

function X = fourierCoeff(t, xt, T, N)
    % Compute Fourier series coefficients numerically
    w0 = 2*pi/T;
    X = zeros(1, 2*N+1);

    for k = -N:N
        integrand = xt .* exp(-1j*k*w0*t);
        xk = (1/T) * trapz(t, integrand);
        X(k+N+1) = xk;
    end
end

function x_hat = partialFourierSum(Ck, T, time_grid)
    % Partial Fourier reconstruction from coefficients
    N = (length(Ck)-1)/2;
    w0 = 2*pi/T;
    x_hat = zeros(size(time_grid));

    for k = -N:N
        x_hat = x_hat + Ck(k+N+1) * exp(1j*k*w0*time_grid);
    end

    x_hat = real(x_hat); 
end
