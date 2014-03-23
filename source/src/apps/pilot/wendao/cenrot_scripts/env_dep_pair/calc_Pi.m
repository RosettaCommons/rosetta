function [ P_i ] = calc_Pi(aa_surf, aa_core, ai, w)
wsurfN = w*aa_surf;
wcoreN = (1-w)*aa_core;
w_N = wsurfN + wcoreN;
sum_r_aa = sum(w_N, 2);
P_i = w_N(ai,:) ./ sum_r_aa(ai);
