syms t C1 C2

a = -2*t + 5;

v = int(a, t) + C1;

v = subs(v, C1, solve(subs(v, t, 0) == 2, C1));

x = int(v, t) + C2;

x = subs(x, C2, solve(subs(x, t, 0) == 5, C2));

a_func = matlabFunction(a);
v_func = matlabFunction(v);
x_func = matlabFunction(x);

t_vals = linspace(0, 5, 500);  

a_vals = a_func(t_vals);
v_vals = v_func(t_vals);
x_vals = x_func(t_vals);

figure;
plot(t_vals, x_vals, 'r-', 'LineWidth', 2); 
hold on;
plot(t_vals, v_vals, 'g', 'LineWidth', 2);
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
hold off;

title('Position, Velocity, and Acceleration vs Time');
xlabel('Time (s)');
ylabel('Magnitude');
legend('Position x(t)', 'Velocity v(t)', 'Acceleration a(t)', 'Location', 'best');
grid on;
