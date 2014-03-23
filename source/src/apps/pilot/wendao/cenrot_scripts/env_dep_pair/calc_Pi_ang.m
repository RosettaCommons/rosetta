function [ P_i ] = calc_Pi_ang(aa_surf, aa_core, ai, w)
wsurfN = w*squeeze(aa_surf(ai,:,:));
wcoreN = (1-w)*squeeze(aa_core(ai,:,:));
w_N = wsurfN + wcoreN;
sum_r_aa = sum(w_N, 2);
P_i = w_N ./ (sum_r_aa*ones(1,16));
