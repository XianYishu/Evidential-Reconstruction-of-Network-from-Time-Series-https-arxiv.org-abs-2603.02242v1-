%% Evidential Network Reconstruction - Main Script
% Based on Dempster-Shafer evidence theory
clear; clc;

% Load ground truth network (Karate club)
load("network-Karate_adj_matrix.mat")
G_graph = Karate_adj_matrix;          % true adjacency matrix
n = length(G_graph);                  % number of nodes

% SIS simulation parameters
Time_max = 2000;                      % length of each time series
num_sources = 5;                      % number of independent infection sources
beta = 0.43;                           % infection rate
gamma = 0.9;                         % recovery rate

% Preallocate cell arrays for BPAs from each source
m_T_cell = cell(1, num_sources);      % belief for connection (m{T})
m_F_cell = cell(1, num_sources);      % belief for disconnection (m{F})
m_TF_cell = cell(1, num_sources);     % belief for uncertainty (m{T,F})
Time_cell = cell(1, num_sources);     % raw time series (for validation)

parfor i = 1:num_sources
    [m_T_cell{i}, m_F_cell{i}, m_TF_cell{i}, Time_cell{i}] = ...
        sis_link_new(G_graph, beta, gamma, Time_max);
end

% Check if any simulation failed (empty output)
hasEmptyMatrix = any(cellfun(@isempty, Time_cell));
if hasEmptyMatrix
    disp('Infection rate too low, please set a higher infection rate');
else
    disp('All simulations completed successfully');
end

%% Multi-source fusion (second-level fusion) using Dempster's rule
% Start with the first source
m_T_final = m_T_cell{1};
m_F_final = m_F_cell{1};
m_TF_final = m_TF_cell{1};

for i = 2:num_sources
    m_T_next = m_T_cell{i};
    m_F_next = m_F_cell{i};
    m_TF_next = m_TF_cell{i};
    
    % Conflict factor K = sum_{X∩Y=∅} m1(X)*m2(Y)
    K = m_T_final .* m_F_next + m_T_next .* m_F_final;
    
    % Dempster's combination rule (Eq. 11-12)
    m_T_final = (m_T_final .* m_T_next + m_T_final .* m_TF_next + m_T_next .* m_TF_final) ./ (1 - K);
    m_F_final = (m_F_final .* m_F_next + m_F_final .* m_TF_next + m_F_next .* m_TF_final) ./ (1 - K);
    m_TF_final = (m_TF_final .* m_TF_next) ./ (1 - K);
    
    % Normalization (ensure sum = 1 for each node pair)
    total = m_T_final + m_F_final + m_TF_final;
    m_T_final = m_T_final ./ total;
    m_F_final = m_F_final ./ total;
    m_TF_final = m_TF_final ./ total;
end

%% Decision Rule Based on Minimum Robustness (DR-MR)
% Find the lowest threshold that yields a connected component covering at least 95% of nodes
epsilon = 1e-8;
r_temp = m_T_final >= epsilon;
G_temp = graph(r_temp);
bins_temp = conncomp(G_temp);
threshold_min = 1;
target_ratio = 0.95;   % require largest component >= 95% of nodes

if sum(bins_temp == 1)/n < target_ratio
    disp('Time series insufficient to reconstruct a fully connected network at default threshold');
else
    threshold_min = 1;
    r_curr = m_T_final >= threshold_min;
    G_curr = graph(r_curr);
    bins_curr = conncomp(G_curr);
    
    for order = 1:8
        step = 10^(-order);
        while sum(bins_curr == 1)/n < target_ratio
            threshold_min = threshold_min - step;
            r_curr = m_T_final >= threshold_min;
            G_curr = graph(r_curr);
            bins_curr = conncomp(G_curr);
            if threshold_min < 0
                threshold_min = 0;
                break;
            end
        end
        if sum(bins_curr == 1)/n < target_ratio && threshold_min ~= 0
            break
        end
        threshold_min = threshold_min + step;
        r_curr = m_T_final >= threshold_min;
        G_curr = graph(r_curr);
        bins_curr = conncomp(G_curr);
    end
end

%% Decision Rule Based on Maximum Similarity (DR-MS)
% Find threshold that maximizes Jaccard similarity between predicted infections and real infections
threshold_range = linspace(1e-8, threshold_min, 400);
Jaccard_vals = zeros(length(threshold_range), 1);

parfor idx = 1:length(threshold_range)
    rho = threshold_range(idx);
    A_recon = m_T_final >= rho;
    sim_vals = zeros(1, 5);
    for j = 1:5
        sim_vals(j) = compute_infection_jaccard(A_recon, Time_cell{randi(num_sources)}, 1000);
    end
    Jaccard_vals(idx) = mean(sim_vals);
end

% Select threshold that gives maximum Jaccard similarity
optimal_rho = threshold_range(Jaccard_vals == max(Jaccard_vals));
A_recon = m_T_final >= optimal_rho;

% Evaluate reconstruction performance (requires true network for benchmark)
[s, r] = similarity(G_graph, A_recon);
disp(['Reconstruction rate (s): ', num2str(s*100), '%']);
disp(['Redundancy rate (r): ', num2str(r*100), '%']);
