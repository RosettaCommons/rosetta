function [ P_ij ] = calc_Pij(aa_aa, ai, aj)

%smooth
pair_counts = aa_aa;

sum_ij = squeeze(sum(pair_counts(ai, aj, :)));
P_ij = squeeze(pair_counts(ai, aj, :)) ./ sum_ij;
