%% Helper Functions

function [m_T, m_F, m_TF, time_series] = sis_link_new(G_graph, beta, gamma, Time_max, seed_nodes)
% Simulate SIS dynamics and compute first-level fused BPA for a single source.
% Input:
%   G_graph  - n x n true adjacency matrix
%   beta     - infection rate
%   gamma    - recovery rate
%   Time_max - maximum time steps
% Output:
%   m_T, m_F, m_TF - fused BPA matrices (belief for connection, disconnection, uncertainty)
%   time_series    - simulated binary time series (Tmax x n)
    n = length(G_graph);
    time_series = zeros(Time_max, n);
    % fixed initial infected nodes
    if ~exist('seed_nodes', 'var') || isempty(seed_nodes)
        seed_nodes = zeros(1, 10);  % 预分配
        for s = 1:10
            seed_nodes(s) = randi(n);
        end
    end
    time_series(1, seed_nodes) = 1;
    
    % Precompute neighbor lists for efficiency
    neighbors = cell(1, n);
    for i = 1:n
        neighbors{i} = find(G_graph(i, :) == 1);
    end
    
    % Repeated simulation until infection survives to maximum time
    attempt = 0;
    while true
        attempt = attempt + 1;
        if attempt > 500
            m_T = []; m_F = []; m_TF = []; time_series = [];
            return;
        end
        
        T = 1;
        state = zeros(1, n);
        state(seed_nodes) = 1;
        time_series(1, :) = state;
        
        while true
            T = T + 1;
            infected_now = find(state == 1);
            
            % Collect all neighbors of currently infected nodes
            neighbor_list = [];
            for i = 1:length(infected_now)
                neighbor_list = [neighbor_list, neighbors{infected_now(i)}];
            end
            neighbor_list = unique(neighbor_list);
            neighbor_list(randperm(length(neighbor_list)));  % random permutation (not used further)
            
            % Infection attempts
            for i = 1:length(neighbor_list)
                if rand <= beta
                    state(neighbor_list(i)) = 1;
                end
            end
            
            % Recovery attempts
            infected_prev = find(time_series(T-1, :) == 1);
            for i = 1:length(infected_prev)
                if rand <= gamma
                    state(infected_prev(i)) = 0;
                end
            end
            
            time_series(T, :) = state;
            
            % Stopping conditions
            if all(state == 0)
                break;
            elseif T == Time_max
                break;
            end
        end
        
        if T == Time_max
            break;
        end
    end
    
    %% Build association matrix M_A and non-association matrix M_N
    M_A = zeros(n, n);
    M_N = zeros(n, n);
    
    for t = 1:(size(time_series,1)-1)
        infected_t = find(time_series(t,:) == 1);
        susceptible_t = find(time_series(t,:) == 0);
        susceptible_t1 = find(time_series(t+1,:) == 0);
        infected_t1 = find(time_series(t+1,:) == 1);
        
        S00 = intersect(susceptible_t, susceptible_t1);   % 0->0
        S01 = intersect(susceptible_t, infected_t1);      % 0->1
        S10 = intersect(infected_t, susceptible_t1);      % 1->0
        S11 = intersect(infected_t, infected_t1);         % 1->1
        
        % Update M_N (non-association) - Eq.(4)
        w = [S11, S10, S01];
        M_N(S00, w) = M_N(S00, w) + 1;
        M_N(w, S00) = M_N(w, S00) + 1;
        M_N(S01, S01) = M_N(S01, S01) + 1;
        
        % Update M_A (association) - Eq.(3)
        M_A(S10, S01) = M_A(S10, S01) + 1;
        M_A(S01, S10) = M_A(S01, S10) + 1;
        M_A(S11, S01) = M_A(S11, S01) + 1;
        M_A(S01, S11) = M_A(S01, S11) + 1;
    end
    
    %% Generate BPA from M_N (non-association) - Eqs.(8)-(10)
    M_N(logical(eye(size(M_N)))) = NaN;   % ignore self-loops
    max_N = max(M_N, [], 'all');
    min_N = min(M_N, [], 'all');
    mu_N = (min_N + max_N) / 2;
    
    m_N_T = (mu_N - M_N) ./ (max_N - mu_N);   % belief for connection (from non-association)
    m_N_T(m_N_T < 0) = 0;
    m_N_F = (M_N - mu_N) ./ (mu_N - min_N);   % belief for disconnection
    m_N_F(m_N_F < 0) = 0;
    m_N_TF = zeros(n);
    for i = 1:n
        for j = 1:n
            if M_N(i,j) >= mu_N
                m_N_TF(i,j) = (max_N - M_N(i,j)) / (max_N - mu_N);
            else
                m_N_TF(i,j) = (M_N(i,j) - min_N) / (mu_N - min_N);
            end
        end
    end
    
    %% Generate BPA from M_A (association) - Eqs.(5)-(7)
    M_A(logical(eye(size(M_A)))) = NaN;
    max_A = max(M_A, [], 'all');
    min_A = min(M_A, [], 'all');
    mu_A = (min_A + max_A) / 2;
    
    m_A_T = (M_A - mu_A) ./ (mu_A - min_A);   % belief for connection (from association)
    m_A_T(m_A_T < 0) = 0;
    m_A_F = (mu_A - M_A) ./ (max_A - mu_A);   % belief for disconnection
    m_A_F(m_A_F < 0) = 0;
    m_A_TF = zeros(n);
    for i = 1:n
        for j = 1:n
            if M_A(i,j) >= mu_A
                m_A_TF(i,j) = (max_A - M_A(i,j)) / (max_A - mu_A);
            else
                m_A_TF(i,j) = (M_A(i,j) - min_A) / (mu_A - min_A);
            end
        end
    end
    
    %% First-level fusion: combine M_A and M_N BPAs using Dempster's rule
    K = m_A_T .* m_N_F + m_N_T .* m_A_F;
    m_T = (m_A_T .* m_N_T + m_A_T .* m_N_TF + m_N_T .* m_A_TF) ./ (1 - K);
    m_F = (m_A_F .* m_N_F + m_A_F .* m_N_TF + m_N_F .* m_A_TF) ./ (1 - K);
    m_TF = (m_A_TF .* m_N_TF) ./ (1 - K);
end

