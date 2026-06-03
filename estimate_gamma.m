gamma_estimated = estimate_gamma_from_time_series(Time_cell);
disp(['Estimated recovery rate gamma = ', num2str(gamma_estimated)]);

function gamma_est = estimate_gamma_from_time_series(Ts_cell)
% Estimate recovery rate gamma from multi-source SIS time series.
% Input:
%   Ts_cell - cell array, each cell is a Tmax x n binary matrix (0=S,1=I)
% Output:
%   gamma_est - estimated recovery rate (scalar)
%
% Based on Eq. in Section III (paper page 10):
%   r_t = 1 - |V_inf(t+1) ∩ V_inf(t)| / |V_inf(t)|
%   E[r_t] = gamma
%   gamma_est = average of r_t over all time steps and all sources.

    n_sources = length(Ts_cell);
    total_r = [];
    
    for s = 1:n_sources
        Ts = Ts_cell{s};
        [Tmax, n] = size(Ts);
        if Tmax < 2
            continue;
        end
        for t = 1:Tmax-1
            infected_t = find(Ts(t,:) == 1);
            infected_t1 = find(Ts(t+1,:) == 1);
            n_inf_t = length(infected_t);
            if n_inf_t == 0
                continue;   % skip if no infected nodes at time t
            end
            % Number of nodes that remain infected from t to t+1
            n_still_infected = sum(ismember(infected_t, infected_t1));
            r_t = 1 - n_still_infected / n_inf_t;
            total_r = [total_r, r_t];
        end
    end
    
    if isempty(total_r)
        gamma_est = NaN;
        warning('No valid time steps found for gamma estimation');
    else
        gamma_est = mean(total_r);
    end
end