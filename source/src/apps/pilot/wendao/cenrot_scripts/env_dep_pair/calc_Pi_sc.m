function [ P_i ] = calc_Pi_sc(ss, sc, cs, cc, wi, wj)
wss = wi*wj*squeeze(ss);
wsc = wi*(1-wj)*squeeze(sc);
wcs = (1-wi)*wj*squeeze(cs);
wcc = (1-wi)*(1-wj)*squeeze(cc);
w_N = squeeze(wss + wsc + wcs + wcc);
sum_r_aa = sum(w_N);
P_i = w_N ./ sum_r_aa;

