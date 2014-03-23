function [ P_i ] = calc_NtoP(N)
N = squeeze(N);
sum_r_aa = sum(N, 2);
P_i = N ./ (sum_r_aa*ones(1,16));
P_i(1:8,:)=0.0625;

