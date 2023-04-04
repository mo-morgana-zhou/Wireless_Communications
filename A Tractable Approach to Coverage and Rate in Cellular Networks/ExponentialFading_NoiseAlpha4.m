%% Exponential Fading, Noise, 𝛼 = 4
% Simulation parameters
service_provider_index = 1;
tx_power = 1; % transmitted power
mu = 1;
sigma2 = 1e-7; % sigma2: noise power
path_loss_exponent = 4; % path loss exponents

%% Genrate the real world scenario
L = 1000; % length of the square region
Area = L^2;
lambda = 1e-4; % intensity
numbPoints = 100000;

%% Calculate the SNR and SINR
SINR_values = zeros(numbPoints, 1);
SimulationTimes = 10000;

for i = 1 : SimulationTimes
    % Generate base station locations
    base_station_locations = L * rand(numbPoints, 2);
    
    % Generate mobile user locations
    mobile_user_locations = [0, 0];
    
    % Associate each mobile user with the closest base station
    user_bs_indices = dsearchn(base_station_locations, mobile_user_locations);
    
    % Calculate the distance from the mobile user to its associated base station
    distance = norm(base_station_locations(user_bs_indices,:) - mobile_user_locations);
    
    % Calculate the received power at the mobile user from its associated base station
    signal_power = tx_power * distance^(-path_loss_exponent);
    
    % Calculate the interference power at the mobile user from the interfering base stations
    bs_distances = sqrt(sum((base_station_locations - base_station_locations(user_bs_indices,:)).^2, 2)); % calculate distances between base stations
    bs_distances(user_bs_indices) = inf; % set distance to associated base station to infinity
    interference_power = sum(tx_power * bs_distances.^(-path_loss_exponent)); % sum interference power over interfering base stations
    
    % Calculate the SINR at the mobile user
    SINR_values(i) = signal_power / (interference_power + sigma2);        
end

%% Function to calculate the coverage probability in real-world scenario
coverage_prob = @(T) sum(SINR_values>T) / SimulationTimes;

% Threshold values
T = logspace(-1, 1, 10000);

% Calculate coverage probabilities
coverage_prob_values = arrayfun(@(x) coverage_prob(x), T);

%% Calculate p_c using provided function
SNR = 1/(mu * sigma2);

% Function to calculate: rho(T,4) = √𝑇(𝜋/2−arctan(1/√T)).
rho = @(T) T^(1/2) * (pi/2 - atan(1/T^(1/2)));

% Function to calculate: kappa(T) = 1 + rho(T,4)
kappa = @(T) 1 + rho(T);

% Function to calculate the coverage probability
p_c = @(T, lambda, SNR) (pi^(3/2)*lambda)/(sqrt(T/SNR)) * exp((lambda*pi*kappa(T))^2 / (4*T/SNR)) ...
    * qfunc((lambda*pi*kappa(T)) / sqrt(2*T/SNR));

% Calculate coverage probabilities
p_c_values = arrayfun(@(x) p_c(x, lambda, SNR), T); % To keep the dimensions consistent with T (array)

%% Plot the function
% Plot the coverage probability from the real world simulation
semilogx(T, coverage_prob_values);
hold on;

% Plot the function to calculate coverage probability
semilogx(T, p_c_values);

xlabel('Threshold');
ylabel('Coverage Probability');
title('Coverage Probability vs Threshold');
legend('Real World Simulation', 'Function Calculation');


