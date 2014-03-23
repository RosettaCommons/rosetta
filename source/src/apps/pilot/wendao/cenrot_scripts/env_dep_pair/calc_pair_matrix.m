function [ matrix ] = calc_pair_matrix(pair)
na = sum(pair);
all = sum(na);
p_a = na / all;
papa = p_a' * p_a;
paa = pair / all;
matrix = paa ./ papa; 
