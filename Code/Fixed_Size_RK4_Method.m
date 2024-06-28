clear all;
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta 4 stage's method for the IVP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0 = 3;
Gamma = 1/14;
Sigma = 1/7;
Omega = 1/365;
Mue = 1/(76*365);
Alpha = 0;
Beta = R0 * Gamma;

tend = 2000; % final time

% Initial condition of ODE's
S0 = 0.999;
E0 = 0.001;
I0 = 0.000;
R0 = 0.000;

time_steps = [0.01, 0.1, 1, 2, 5, 10]; % different time steps

% Loop over different time steps
for j = 1:length(time_steps)
    h = time_steps(j);
    
    % Initialize time and solution vectors
    time = 0:h:tend;
    y = zeros(4, length(time));
    y(:,1) = [S0; E0; I0; R0]; % Initial values
    
    tic; % Start timing fixed step size method
    for i = 1:length(time)-1
        k1 = ODESystem(time(i), y(:,i), Alpha, Beta, Gamma, Sigma, Mue, Omega);
        k2 = ODESystem(time(i)+h/2, y(:,i)+h*k1/2, Alpha, Beta, Gamma, Sigma, Mue, Omega);
        k3 = ODESystem(time(i)+h/2, y(:,i)+h*k2/2, Alpha, Beta, Gamma, Sigma, Mue, Omega);
        k4 = ODESystem(time(i)+h, y(:,i)+h*k3, Alpha, Beta, Gamma, Sigma, Mue, Omega);
        y(:,i+1) = y(:,i) + (1/6)*(k1+2*k2+2*k3+k4)*h;
    end
    fixed_time = toc; % End timing fixed step size method
    
    S = y(1,:);
    E = y(2,:);
    I = y(3,:);
    R = y(4,:);
    
    disp(['Fixed Step Size Method - Execution Time for h = ', num2str(h), ': ', num2str(fixed_time), ' seconds']);
    disp(['Final values for h = ', num2str(h), ': S = ', num2str(S(end)), ', E = ', num2str(E(end)), ', I = ', num2str(I(end)), ', R = ', num2str(R(end))]);
    
    figure;
    plot(time, S, 'k', 'DisplayName', 'Susceptible');
    hold on;
    plot(time, E, 'g', 'DisplayName', 'Exposed');
    plot(time, I, 'r', 'DisplayName', 'Infected');
    plot(time, R, 'b', 'DisplayName', 'Recovered');
    xlabel('Time [days]');
    ylabel('Population');
    legend('show');
    title(['SEIR Model Simulation using Fixed Step Size Method (h = ', num2str(h), ')']);
    hold off;
end

% Nested function for ODE system
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
