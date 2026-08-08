%% ============================================================
%  Inverted Pendulum Stabilization via LQR Optimal Control
%  ============================================================
%  Model  : Cart-Pole System (Inverted Pendulum on a Cart)
%  Input  : Horizontal Force Applied to the Cart (F)
%  Outputs: Cart Position (x) and Pendulum Angle (theta)
%  ============================================================

clear; clc; close all;

%% ============================================================
%  1. Physical Parameters of the System
%  ============================================================

% -- Cart Parameters --
M = 1.0;          % Mass of the cart [kg]

% -- Pendulum Parameters --
m = 0.3;          % Mass of the pendulum bob [kg]
l = 0.5;          % Half-length of the pendulum (pivot to CoM) [m]
g = 9.81;         % Gravitational acceleration [m/s^2]

% -- Friction Coefficients --
b_cart = 0.1;     % Viscous friction coefficient of the cart [N.s/m]
b_pend = 0.01;    % Viscous friction coefficient at the pivot [N.m.s/rad]

% -- Moment of Inertia --
I = m * l^2 / 3;  % Moment of inertia of the pendulum about the pivot [kg.m^2]

fprintf('=== Inverted Pendulum System Parameters ===\n');
fprintf('Cart Mass:          M = %.2f kg\n', M);
fprintf('Pendulum Mass:      m = %.2f kg\n', m);
fprintf('Pendulum Half-Length: l = %.2f m\n', l);
fprintf('Gravity:            g = %.2f m/s^2\n', g);
fprintf('Moment of Inertia:  I = %.4f kg.m^2\n', I);
fprintf('==========================================\n\n');

%% ============================================================
%  2. State-Space Modeling (Linearized about the Upright Equilibrium)
%  ============================================================
%  State vector: x = [x; x_dot; theta; theta_dot]
%    x         : Cart position [m]
%    x_dot     : Cart velocity [m/s]
%    theta     : Pendulum angle from the vertical [rad]
%    theta_dot : Angular velocity of the pendulum [rad/s]
%
%  Linearized equations of motion about theta = 0 (upright):
%    (M+m)*x_dd + m*l*theta_dd = F - b_cart*x_d
%    m*l*x_dd + (I+m*l^2)*theta_dd = m*g*l*theta - b_pend*theta_d
%
%  Solving for the accelerations via mass matrix inversion:
%    x_dd     = [(I+ml^2)(F-b_c*x_d) - ml(mgl*theta - b_p*theta_d)] / D
%    theta_dd = [(M+m)(mgl*theta - b_p*theta_d) - ml(F-b_c*x_d)]   / D
%
%  Resulting gravity-coupling signs:
%    A(2,3) = -m^2*g*l^2 / D   (negative)
%    A(4,3) = +(M+m)*m*g*l / D  (positive => inherent instability)

% System matrices
Mt = M + m;
D = Mt * (I + m*l^2) - (m*l)^2;

% State matrix A
A = [0,             1,                    0,             0;             ...
     0,  -(I+m*l^2)*b_cart/D,  -(m^2*g*l^2)/D,   m*l*b_pend/D;  ...
     0,             0,                    0,             1;             ...
     0,  (m*l*b_cart)/D,        (Mt*m*g*l)/D,   -(Mt*b_pend)/D];

% Input matrix B
B = [0;                                    ...
     (I + m*l^2)/D;                       ...
     0;                                    ...
     -(m*l)/D];

% Output matrix C (measured outputs: cart position and pendulum angle)
C = [1, 0, 0, 0;                          ...
     0, 0, 1, 0];

% Feedthrough matrix
D_mat = [0; 0];

% Construct the state-space model
sys = ss(A, B, C, D_mat);

fprintf('=== State-Space Model ===\n');
fprintf('System order: %d\n', size(A, 1));
fprintf('Outputs: x (cart position), theta (pendulum angle)\n\n');

% Controllability and observability analysis
Co = ctrb(sys);
Ob = obsv(sys);

fprintf('Controllability matrix rank: %d (expected: %d)\n', rank(Co), size(A,1));
fprintf('Observability matrix rank:  %d (expected: %d)\n\n', rank(Ob), size(A,1));

if rank(Co) == size(A,1)
    fprintf('>>> System is fully controllable.\n');
else
    fprintf('>>> WARNING: System is NOT controllable!\n');
end

if rank(Ob) == size(A,1)
    fprintf('>>> System is fully observable.\n\n');
else
    fprintf('>>> WARNING: System is NOT observable!\n\n');
end

%% ============================================================
%  3. Open-Loop Stability Analysis
%  ============================================================

eig_open = eig(A);
fprintf('=== Open-Loop Poles (Eigenvalues of A) ===\n');
for i = 1:length(eig_open)
    fprintf('  Pole %d: %.4f %+.4fi\n', i, real(eig_open(i)), imag(eig_open(i)));
end

unstable = any(real(eig_open) > 0);
if unstable
    fprintf('>>> The open-loop system is UNSTABLE (at least one RHP pole exists).\n\n');
else
    fprintf('>>> The open-loop system is stable.\n\n');
end

%% ============================================================
%  4. LQR Controller Design
%  ============================================================
%  Cost function: J = integral_0^inf (x'*Q*x + u'*R*u) dt
%  Control law:   u = -K*x

% Weighting matrices Q and R
% Q penalizes state deviations (larger => tighter regulation on that state)
% R penalizes control effort   (larger => more energy-efficient, slower response)

Q = [100,   0,   0,   0;    % Cart position
       0,   1,   0,   0;    % Cart velocity
       0,   0, 500,   0;    % Pendulum angle (highest priority)
       0,   0,   0,  50];   % Pendulum angular velocity

R = 1;                        % Control effort penalty

% Solve the algebraic Riccati equation
[K, S, eigs_cl] = lqr(A, B, Q, R);

fprintf('=== LQR Controller Results ===\n');
fprintf('State-feedback gain vector K:\n');
fprintf('  K = [%.4f, %.4f, %.4f, %.4f]\n', K(1), K(2), K(3), K(4));
fprintf('\nClosed-loop poles:\n');
for i = 1:length(eigs_cl)
    fprintf('  Pole %d: %.4f %+.4fi\n', i, real(eigs_cl(i)), imag(eigs_cl(i)));
end

if all(real(eigs_cl) < 0)
    fprintf('>>> The closed-loop system is STABLE under LQR control.\n\n');
else
    fprintf('>>> WARNING: The closed-loop system is UNSTABLE!\n\n');
end

%% ============================================================
%  5. Time-Domain Simulation
%  ============================================================

% Initial condition: pendulum slightly displaced from the upright position
x0 = [0;           % Initial cart position [m]
      0;           % Initial cart velocity [m/s]
      0.2;         % Initial pendulum angle [rad] (~11.5 deg)
      0];          % Initial angular velocity [rad/s]

% Closed-loop system: u = -K*x
A_cl = A - B * K;
sys_cl = ss(A_cl, B, C, D_mat);

% Simulation time vector
t_end = 10;
dt = 0.01;
t = 0:dt:t_end;

% Open-loop simulation (for comparison)
sys_ol = ss(A, B, C, D_mat);
[y_ol, t_ol, x_ol] = lsim(sys_ol, zeros(length(t), 1), t, x0);

% Closed-loop LQR simulation
[y_lqr, t_lqr, x_lqr] = lsim(sys_cl, zeros(length(t), 1), t, x0);

% Compute the control signal: u = -K*x
u_lqr = -K * x_lqr';

fprintf('=== Simulation Results ===\n');
fprintf('Simulation time: 0 to %d s\n', t_end);
fprintf('Initial conditions: x=%.1f m, x_dot=%.1f m/s, theta=%.2f rad, theta_dot=%.1f rad/s\n', ...
    x0(1), x0(2), x0(3), x0(4));
fprintf('Steady-state cart position error: %.6f m\n', y_lqr(end, 1));
fprintf('Steady-state pendulum angle error: %.6f rad\n\n', y_lqr(end, 2));

%% ============================================================
%  6. Simulation Plots
%  ============================================================

figure('Name', 'Inverted Pendulum - Simulation Results', ...
       'Position', [100, 100, 1200, 800]);

% --- Subplot 1: Cart position ---
subplot(3, 2, 1);
plot(t_ol, y_ol(:, 1), 'r--', 'LineWidth', 1.5); hold on;
plot(t_lqr, y_lqr(:, 1), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cart Position x (m)');
title('Cart Position');
legend('Open-Loop', 'LQR Controlled');
grid on;

% --- Subplot 2: Pendulum angle ---
subplot(3, 2, 2);
plot(t_ol, y_ol(:, 2) * 180/pi, 'r--', 'LineWidth', 1.5); hold on;
plot(t_lqr, y_lqr(:, 2) * 180/pi, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Pendulum Angle (deg)');
title('Pendulum Angle');
legend('Open-Loop', 'LQR Controlled');
grid on;

% --- Subplot 3: Cart velocity ---
subplot(3, 2, 3);
plot(t_lqr, x_lqr(:, 2), 'g-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Cart Velocity (m/s)');
title('Cart Velocity (LQR)');
grid on;

% --- Subplot 4: Angular velocity ---
subplot(3, 2, 4);
plot(t_lqr, x_lqr(:, 4) * 180/pi, 'm-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Velocity (deg/s)');
title('Pendulum Angular Velocity (LQR)');
grid on;

% --- Subplot 5: Control signal ---
subplot(3, 2, 5);
plot(t_lqr, u_lqr, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Force F (N)');
title('Control Signal (Force Applied to Cart)');
grid on;

% --- Subplot 6: Pole-zero map ---
subplot(3, 2, 6);
plot(real(eig_open), imag(eig_open), 'rx', 'MarkerSize', 12, 'LineWidth', 2); hold on;
plot(real(eigs_cl), imag(eigs_cl), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
xline(0, 'k--', 'LineWidth', 1);
yline(0, 'k--', 'LineWidth', 1);
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Pole Map: Open-Loop vs. Closed-Loop');
legend('Open-Loop', 'Closed-Loop (LQR)');
axis equal;
grid on;

sgtitle('Inverted Pendulum Stabilization via LQR', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================
%  7. Animation of the Inverted Pendulum
%  ============================================================

figure('Name', 'Inverted Pendulum Animation', 'Position', [150, 150, 800, 400]);

% Plot range
x_min = min(x_lqr(:,1)) - 1;
x_max = max(x_lqr(:,1)) + 1;

% Initialize axes
axis([x_min, x_max, -0.4, 2*l + 0.3]);
axis equal;
xlabel('Position (m)');
ylabel('Height (m)');
grid on;
hold on;

% Ground
plot([x_min, x_max], [-0.2, -0.2], 'k-', 'LineWidth', 2);

% Wheel templates (fixed geometry, only position changes)
theta_wheel = 0:0.1:2*pi;
xL0 = -0.15 + 0.05*cos(theta_wheel);
yL0 = -0.15 - 0.05 + 0.05*sin(theta_wheel);
xR0 =  0.15 + 0.05*cos(theta_wheel);
yR0 = yL0;

% Create graphical objects (handles for efficient update)
h_cart = rectangle('Position', [0-0.3, -0.15, 0.6, 0.15], ...
    'FaceColor', [0.2, 0.4, 0.8], 'EdgeColor', 'k', 'LineWidth', 2);
h_wL = plot(xL0, yL0, 'k-', 'LineWidth', 1.5);
h_wR = plot(xR0, yR0, 'k-', 'LineWidth', 1.5);
h_rod = plot([0, 0], [0, 2*l], 'r-', 'LineWidth', 3);
h_axl = plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
h_bob = plot(0, 2*l, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', [0.8, 0.2, 0.2]);
h_title = title('', 'FontSize', 12);

% Frame skip factor for faster rendering
skip = 3;
for i = 1:skip:length(t)
    x_cart = x_lqr(i, 1);
    theta  = x_lqr(i, 3);

    x_pend = x_cart + 2*l * sin(theta);
    y_pend = 2*l * cos(theta);

    % Update graphical handles (avoids costly clf + redraw)
    set(h_cart, 'Position', [x_cart-0.3, -0.15, 0.6, 0.15]);
    set(h_wL, 'XData', x_cart + xL0, 'YData', yL0);
    set(h_wR, 'XData', x_cart + xR0, 'YData', yR0);
    set(h_rod, 'XData', [x_cart, x_pend], 'YData', [0, y_pend]);
    set(h_axl, 'XData', x_cart, 'YData', 0);
    set(h_bob, 'XData', x_pend, 'YData', y_pend);
    set(h_title, 'String', sprintf('Time: %.2f s  |  Angle: %.2f deg  |  Force: %.2f N', ...
          t(i), theta*180/pi, u_lqr(i)));

    drawnow;
end
hold off;
fprintf('Animation completed.\n');

%% ============================================================
%  8. Sensitivity Analysis - Effect of Pendulum Mass Variation
%  ============================================================

fprintf('\n=== Sensitivity Analysis ===\n');

mass_values = [0.1, 0.3, 0.5, 1.0, 2.0];
fprintf('\nEffect of pendulum mass on closed-loop poles:\n');
fprintf('%-12s %-28s %-28s\n', 'Mass (kg)', 'Pole 1', 'Pole 2');
fprintf('%s\n', repmat('-', 1, 72));

figure('Name', 'Sensitivity Analysis', 'Position', [200, 200, 1000, 400]);

for idx = 1:length(mass_values)
    m_test = mass_values(idx);
    I_test = m_test * l^2 / 3;
    Mt_test = M + m_test;
    D_test = Mt_test * (I_test + m_test*l^2) - (m_test*l)^2;
    
    A_test = [0,             1,                          0,                    0;  ...
               0,  -(I_test+m_test*l^2)*b_cart/D_test,  -(m_test^2*g*l^2)/D_test,   m_test*l*b_pend/D_test;  ...
               0,             0,                          0,                    1;  ...
               0,  (m_test*l*b_cart)/D_test,         (Mt_test*m_test*g*l)/D_test,   -(Mt_test*b_pend)/D_test];
    
    B_test = [0; (I_test + m_test*l^2)/D_test; 0; -(m_test*l)/D_test];
    
    [K_test, ~, eigs_test] = lqr(A_test, B_test, Q, R);
    
    fprintf('%-12.2f ', m_test);
    for j = 1:min(2, length(eigs_test))
        fprintf('%-28s ', sprintf('%.3f %+.3fi', real(eigs_test(j)), imag(eigs_test(j))));
    end
    fprintf('\n');
    
    % Simulate with the modified mass
    x0_test = [0; 0; 0.2; 0];
    A_cl_test = A_test - B_test * K_test;
    sys_cl_test = ss(A_cl_test, B_test, C, D_mat);
    [y_test, ~, ~] = lsim(sys_cl_test, zeros(length(t), 1), t, x0_test);
    
    subplot(1, 2, 1);
    plot(t, y_test(:, 2) * 180/pi, 'LineWidth', 1.5); hold on;
    
    subplot(1, 2, 2);
    plot(t, y_test(:, 1), 'LineWidth', 1.5); hold on;
end

subplot(1, 2, 1);
xlabel('Time (s)'); ylabel('Pendulum Angle (deg)');
title('Effect of Pendulum Mass on Angle');
legend(arrayfun(@(x) sprintf('m=%.1f kg', x), mass_values, 'UniformOutput', false));
grid on;

subplot(1, 2, 2);
xlabel('Time (s)'); ylabel('Cart Position (m)');
title('Effect of Pendulum Mass on Cart Position');
legend(arrayfun(@(x) sprintf('m=%.1f kg', x), mass_values, 'UniformOutput', false));
grid on;

sgtitle('Sensitivity Analysis with Respect to Pendulum Mass', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================
%  9. Comparative Study: LQR vs. PID Controller
%  ============================================================

fprintf('\n=== Comparative Study: LQR vs. PID ===\n');

% Derive the pendulum angle transfer function from the state-space model
% G_theta(s) = C_theta * (sI - A)^{-1} * B
C_theta = [0, 0, 1, 0];
sys_theta = ss(A, B, C_theta, 0);
G_theta = tf(sys_theta);

% PID controller gains (tuned for the angle subsystem)
Kp_pid = -200;
Ki_pid = -50;
Kd_pid = -40;
C_pid = pid(Kp_pid, Ki_pid, Kd_pid);

% Closed-loop PID system
sys_pid_cl = feedback(C_pid * G_theta, 1);

% Simulate PID step response
[y_pid, t_pid] = lsim(sys_pid_cl, 0.2*ones(length(t), 1), t);

% Compute the PID control signal
e_pid = 0.2*ones(length(t), 1) - y_pid;
u_pid = zeros(size(t));
for i = 2:length(t)
    u_pid(i) = Kp_pid*e_pid(i) + Ki_pid*trapz(t(1:i), e_pid(1:i)) + Kd_pid*(e_pid(i)-e_pid(i-1))/dt;
end

figure('Name', 'LQR vs. PID Comparison', 'Position', [250, 250, 800, 400]);

subplot(1, 2, 1);
plot(t, y_lqr(:, 2) * 180/pi, 'b-', 'LineWidth', 2); hold on;
plot(t_pid, y_pid * 180/pi, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Pendulum Angle (deg)');
title('Angle Response Comparison');
legend('LQR', 'PID');
grid on;

subplot(1, 2, 2);
plot(t, u_lqr, 'b-', 'LineWidth', 2); hold on;
plot(t_pid, u_pid, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Control Force (N)');
title('Control Effort Comparison');
legend('LQR', 'PID');
grid on;

sgtitle('Comparative Study: LQR vs. PID Controller', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================
%  10. Save Results
%  ============================================================

save('inverted_pendulum_results.mat', 'A', 'B', 'C', 'K', 'Q', 'R', ...
     'sys', 'sys_cl', 'x_lqr', 'y_lqr', 'u_lqr', 't', ...
     'eigs_cl', 'eig_open');

fprintf('\n=== Simulation Complete ===\n');
fprintf('Results saved to inverted_pendulum_results.mat\n');
fprintf('Plots and animation have been displayed.\n');
