%% KL Validation: Find the closest pair between original-original and simulated-original KL divergences
% This script computes KL divergences (relative entropy) between:
%   1. All pairs of original multi-source time series (baseline D_KL),
%   2. Each simulated time series (generated from reconstructed network with same seeds) against all original series (hat_D_KL).
% Then selects the closest values from the two sets and evaluates reliability via Eq.(18).
%
% Input variables (from main script):
%   Time_cell               - cell array of original time series (Tmax x n binary)
%   A_recon                 - reconstructed adjacency matrix (n x n)
%   estimated_beta          - estimated infection rate (gamma_estimated)
%   gamma_estimated         - estimated recovery rate (optimal_rho)
%   Time_max                - length of each time series (Tmax)
%
% Output: displays relative deviation delta_KL as "Redundancy rate (r): ..."

num_sources = length(Time_cell);
num_sim = num_sources;   % one simulated series per original source (seed-synchronized)

% 1. Compute all pairwise KL divergences among original series (i<j)
kl_orig = [];
for i = 1:num_sources
    for j = i+1:num_sources
        kl_val = KL(Time_cell{i}, Time_cell{j});
        kl_orig = [kl_orig, kl_val];
    end
end

% 2. For each original source, simulate a time series using the same initial infected seeds,
%    then compute KL between this simulated series and every original series.
kl_sim = [];
for idx = 1:num_sim
    % extract initial infected nodes from the idx-th original series
    seed_nodes = find(Time_cell{idx}(1, :) == 1);
    % simulate with identical seeds
    time_sim = simulate_sis_with_seeds(A_recon, estimated_beta, gamma_estimated, Time_max, seed_nodes);
    % compute KL with all original series
    for j = 1:num_sources
        kl_val = KL(Time_cell{j}, time_sim);
        kl_sim = [kl_sim, kl_val];
    end
end

% 3. Find the pair (one from kl_orig, one from kl_sim) with the smallest absolute difference
min_diff = inf;
best_orig = NaN;
best_sim = NaN;
for i = 1:length(kl_orig)
    for j = 1:length(kl_sim)
        diff = abs(kl_orig(i) - kl_sim(j));
        if diff < min_diff
            min_diff = diff;
            best_orig = kl_orig(i);
            best_sim = kl_sim(j);
        end
    end
end

% 4. Apply the decision rule (Eq.18): delta = |D_KL - hat_D_KL| / D_KL
kk = best_orig;   % baseline D_KL (smallest among original pairs? Actually here it is the closest one, but variable name kept)
k = best_sim;     % simulated hat_D_KL
Delta_k = abs(k - kk) / k;
disp(['DETA_KL: ', num2str(Delta_k*100), '%']);

%% Helper function: SIS simulation with given initial infected seeds
function time_series = simulate_sis_with_seeds(G_graph, beta, gamma, Time_max, seed_nodes)
    n = length(G_graph);
    time_series = zeros(Time_max, n);
    time_series(1, seed_nodes) = 1;
    
    neighbors = cell(1, n);
    for i = 1:n
        neighbors{i} = find(G_graph(i, :) == 1);
    end
    
    state = zeros(1, n);
    state(seed_nodes) = 1;
    time_series(1, :) = state;
    
    for T = 2:Time_max
        infected_now = find(state == 1);
        neighbor_list = [];
        for i = infected_now
            neighbor_list = [neighbor_list, neighbors{i}];
        end
        neighbor_list = unique(neighbor_list);
        
        for i = neighbor_list
            if rand <= beta
                state(i) = 1;
            end
        end
        
        infected_prev = find(time_series(T-1, :) == 1);
        for i = infected_prev
            if rand <= gamma
                state(i) = 0;
            end
        end
        
        time_series(T, :) = state;
        
        if all(state == 0)
            if T < Time_max
                time_series(T+1:end, :) = 0;
            end
            break;
        end
    end
end

%% KL divergence function: uses empirical infection probability distribution (no min-max scaling)
function kl = KL(time1, time2)
    freq1 = sum(time1, 1);
    freq2 = sum(time2, 1);
    prob1 = freq1 / sum(freq1);
    prob2 = freq2 / sum(freq2);
    prob1(prob1 == 0) = eps;
    prob2(prob2 == 0) = eps;
    kl = sum(prob1 .* log(prob1 ./ prob2));
end