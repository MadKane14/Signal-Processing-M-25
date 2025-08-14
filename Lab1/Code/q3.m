function main

nx1 = -10:10; 
x1 = 2*(nx1 == -10) + 2*(nx1 == -10 + 20); % δ[n+10] + δ[n-10]
ny1 = -5:5;
y1 = 3*(ny1 == -5) + 2*(ny1 == 5);         % 3δ[n+5] + 2δ[n-5]
[nz1, z1] = disc_conv(nx1, x1, ny1, y1);

figure;
subplot(3,1,1); stem(nx1, x1); title('x_1[n]'); grid on;
subplot(3,1,2); stem(ny1, y1); title('y_1[n]'); grid on;
subplot(3,1,3); stem(nz1, z1); title('z_1[n] = x_1[n] * y_1[n]'); grid on;

nx2 = -1:1;
x2 = (-1).^nx2;  % (-1)^n
ny2 = -1:1;
y2 = (ny2 == 0) + (ny2 == 1); % δ[n] + δ[n-1]
[nz2, z2] = disc_conv(nx2, x2, ny2, y2);

figure;
subplot(3,1,1); stem(nx2, x2); title('x_2[n]'); grid on;
subplot(3,1,2); stem(ny2, y2); title('y_2[n]'); grid on;
subplot(3,1,3); stem(nz2, z2); title('z_2[n]'); grid on;

nx3 = 0:0; x3 = 4; % constant 4 at n=0
ny3 = -2:1;
y3 = (ny3 == 0) + 2*(ny3 == -1) + (ny3 == -2);
[nz3, z3] = disc_conv(nx3, x3, ny3, y3);

figure;
subplot(3,1,1); stem(nx3, x3); title('x_3[n]'); grid on;
subplot(3,1,2); stem(ny3, y3); title('y_3[n]'); grid on;
subplot(3,1,3); stem(nz3, z3); title('z_3[n]'); grid on;

nx4 = 0:5;
x4 = (nx4 >= 0 & nx4 <= 2); % u[n] - u[n-3]
ny4 = 0:5;
y4 = [0 1 2 ones(1, 3)]; % as per definition
[nz4, z4] = disc_conv(nx4, x4, ny4, y4);

figure;
subplot(3,1,1); stem(nx4, x4); title('x_4[n]'); grid on;
subplot(3,1,2); stem(ny4, y4); title('y_4[n]'); grid on;
subplot(3,1,3); stem(nz4, z4); title('z_4[n]'); grid on;

end

 function [nz, z] = disc_conv(nx, x, ny, y)
    % nx: time indices of x[n]
    % x : samples of x[n]
    % ny: time indices of y[n]
    % y : samples of y[n]

    % Lengths of the signals
    Lx = length(x);
    Ly = length(y);

    % Length of the convolution result
    Lz = Lx + Ly - 1;

    % Initialize output
    z = zeros(1, Lz);

    % Manual convolution using nested loops
    for n = 1:Lz
        for k = 1:Lx
            j = n - k + 1; % index for y
            if j >= 1 && j <= Ly
                z(n) = z(n) + x(k) * y(j);
            end
        end
    end

    % Compute index vector for result
    nz_start = nx(1) + ny(1);
    nz_end   = nx(end) + ny(end);
    nz = nz_start : nz_end;
 end
