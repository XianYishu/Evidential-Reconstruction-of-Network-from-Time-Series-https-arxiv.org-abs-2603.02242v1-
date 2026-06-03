function [reconstruction_rate, redundancy_rate] = similarity(G_true, G_recon)
% Compare reconstructed adjacency with ground truth.
% reconstruction_rate = TP / (TP+FN)
% redundancy_rate = FP / (TP+FN)
    reconstruction_rate = sum(G_true .* G_recon, 'all') / sum(G_true, 'all');
    G_false = 1 - G_true;
    redundancy_rate = sum(G_recon .* G_false, 'all') / sum(G_true, 'all');
end