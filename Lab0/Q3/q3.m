function main
    dt = 0.01;
    t = -5:dt:5;

    h = double(abs(t) <= 1); 

    x1 = double(abs(t) <= 1);           
    x2 = sin(t);                         
    x3 = ((t+1)/2) .* (t >= -1 & t <= 1); 

    [y1, t1] = convolve(t, x1, h, dt);
    [y2, t2] = convolve(t, x2, h, dt);
    [y3, t3] = convolve(t, x3, h, dt);

    figure;
    subplot(2,2,1);
    plot(t,h,'LineWidth',1.5);
    grid on;
    title('Impulse Response h(t)');
    xlabel('t'); ylabel('h(t)');

    subplot(2,2,2);
    plot(t1,y1,'LineWidth',1.5); 
    grid on;
    title('Output for Signal 1.');
    xlabel('t'); ylabel('y(t)');

    subplot(2,2,3);
    plot(t2,y2,'LineWidth',1.5); 
    grid on;
    title('Output for Signal 2.');
    xlabel('t'); ylabel('y(t)');

    subplot(2,2,4);
    plot(t3,y3,'LineWidth',1.5); 
    grid on;
    title('Output for Signal 3.');
    xlabel('t'); ylabel('y(t)');
end

function [y, t_conv] = convolve(t, x, h, dt)
    y = conv(x, h) * dt;  
    
    t_start = t(1) + t(1);             
    t_conv = t_start + (0:length(y)-1) * dt;
end