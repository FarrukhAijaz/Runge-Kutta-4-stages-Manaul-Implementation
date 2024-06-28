clear all;
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta method with adaptive step sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0 = 3;
Gamma = 1/14;
Sigma = 1/7;
Omega = 1/365;
Mue = 1/(76*365);
Alpha = 0;
Beta = R0 * Gamma;

t0 = 0; % initial time
tf = 2000; % final time
h = 0.5; % initial step size

% Tolerances and parameters for adaptive step size control
atol = 1e-6;
rtol = 1e-3;
fac1 = 5;
fac0 = 0.2;
beta = 0.9;
hmin = 0.1;
MAX_ITER = 10000;

% Initial condition of ODE's
S0 = 0.999;
E0 = 0.001;
I0 = 0.000;
R0 = 0.000;

y = [S0;E0;I0;R0];
t = t0;
time = t;
result = y;
step_sizes = h;

ITER = 0;

tic; % Start timing adaptive step size method
while t < tf && ITER < MAX_ITER
    k1 = ODESystem(t, y, Alpha, Beta, Gamma, Sigma, Mue, Omega);
    k2 = ODESystem(t + h/2, y + h*k1/2, Alpha, Beta, Gamma, Sigma, Mue, Omega);
    k3 = ODESystem(t + h/2, y + h*k2/2, Alpha, Beta, Gamma, Sigma, Mue, Omega);
    k4 = ODESystem(t + h, y + h*k3, Alpha, Beta, Gamma, Sigma, Mue, Omega);
    y_new = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
    
    % Error estimation
    eta = (k1 - 2*k2 + 2*k3 - k4) * h/6;
    sigma = norm(eta ./ (atol + rtol * abs(y_new)), 2) / sqrt(length(y_new));
    
    if sigma <= 1
        t = t + h;
        y = y_new;
        time = [time t];
        result = [result y];
        step_sizes = [step_sizes h];
        
        if t >= tf
            break;
        end
    end
    
    % Update step size
    h = h * min(fac1, max(fac0, beta * sigma^(-1/5)));
    if h < hmin
        disp(['Warning: Step size below minimum threshold at iteration ', num2str(ITER), ', t = ', num2str(t)]);
        h = hmin; % Set h to minimum threshold instead of error
    end
    
    ITER = ITER + 1;
end
adaptive_time = toc; % End timing adaptive step size method
S = result(1,:);
E = result(2,:);
I = result(3,:);
R = result(4,:);

disp(['Adaptive Step Size Method - Execution Time: ', num2str(adaptive_time), ' seconds']);
disp(['Final values: S = ', num2str(S(end)), ', E = ', num2str(E(end)), ', I = ', num2str(I(end)), ', R = ', num2str(R(end))]);

figure;
subplot(2,1,1);
plot(time, S, 'k', 'DisplayName', 'Susceptible');
hold on;
plot(time, E, 'g', 'DisplayName', 'Exposed');
plot(time, I, 'r', 'DisplayName', 'Infected');
plot(time, R, 'b', 'DisplayName', 'Recovered');
xlabel('Time [days]');
ylabel('Population');
legend('show');
title('SEIR Model Simulation using Adaptive Runge-Kutta 4th Stage Method');
hold off;

subplot(2,1,2);
plot(time(1:end), step_sizes, 'm', 'DisplayName', 'Step Size');
xlabel('Time [days]');
ylabel('Step Size');
legend('show');
title('Step Sizes over Time');
hold off;

function dydt = ODESystem(time, y, Alpha, Beta, Gamma, Sigma, Mue, Omega)
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    dydt = zeros(4,1);
    dydt(1) = Mue * (S + E + I + R) - (Beta * I * S / (S + E + I + R)) + Omega * R - Mue * S; 
    dydt(2) = (Beta * I * S / (S + E + I + R)) - Sigma * E - Mue * E;
    dydt(3) = Sigma * E - Gamma * I - (Mue + Alpha) * I;
    dydt(4) = Gamma * I - Omega * R - Mue * R;
end
