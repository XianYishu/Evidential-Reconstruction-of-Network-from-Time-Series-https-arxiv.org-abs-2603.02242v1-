function [mean_Jaccard] = compute_infection_jaccard(A_recon, TIME, num_iterations)
% Compute average Jaccard similarity between real new infections (from TIME)
% and predicted new infections (from reconstructed network A_recon).
% This is used in DR-MS to find optimal threshold.
    if nargin < 3
        num_iterations = 2000;
    end
    [Tmax, ~] = size(TIME);
    if ~islogical(TIME)
        TIME_bool = TIME > 0;
    else
        TIME_bool = TIME;
    end
    if ~issparse(A_recon)
        A_recon = sparse(A_recon > 0);
    end
    
    J_vals = zeros(num_iterations, 1);
    max_time = Tmax - 1;
    
    for iter = 1:num_iterations
        t = randi(max_time);
        state_t = TIME_bool(t, :);
        state_t1 = TIME_bool(t+1, :);
        infected_t = find(state_t);
        
        if ~isempty(infected_t)
            % Neighbors of currently infected nodes in reconstructed network
            if isscalar(infected_t)
                neighbor_mask = A_recon(infected_t, :) > 0;
            else
                neighbor_mask = sum(A_recon(infected_t, :), 1) > 0;
            end
            neighbor_mask = neighbor_mask & ~state_t;   % exclude already infected
            
            % Actual new infections (0->1)
            actual_new = state_t1 & ~state_t;
            % Predicted new infections (neighbors that are susceptible)
            predicted_new = neighbor_mask & ~state_t;
            
            inter = sum(actual_new & predicted_new);
            union_val = sum(actual_new | predicted_new);
            if union_val > 0
                J_vals(iter) = inter / union_val;
            end
        end
    end
    
    nonzero_J = J_vals(J_vals > 0);
    if isempty(nonzero_J)
        mean_Jaccard = 0;
    else
        mean_Jaccard = mean(nonzero_J);
    end
end

