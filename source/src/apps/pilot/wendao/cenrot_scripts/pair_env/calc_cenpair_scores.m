function [x, y1, y2]=calc_cenpair_scores(pair_counts,ii,jj,cutoff)
%//counts/2 in load 
pair_counts_symm = pair_counts;
%pair_counts_symm(pair_counts_symm<1) = 0.1; %// avoid log(0)

%smooth data
L=4; %half window
CONV1 = gausswin(2*L+1, 4);
CONV1 = CONV1/sum(CONV1);
for i = 1:20
	for j = i:20
		smoothed = conv(squeeze(pair_counts_symm(i,j,:)), CONV1);
		pair_counts_symm(i,j,:) = smoothed(L+1:end-L);
		if (i!=j) 
			pair_counts_symm(j,i,:) = pair_counts_symm(i,j,:);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P(ai,aj|r)/P(ai|r)/P(aj|r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pair_sums_aan = squeeze(sum(pair_counts_symm)); %20x30
pair_sums_tot = squeeze(sum(pair_sums_aan));    %1x30
totalsum = sum(pair_sums_tot);

%//cal P_ii and P_jj
P_ii = pair_sums_aan(ii, :) ./ pair_sums_tot;
P_jj = pair_sums_aan(jj, :) ./ pair_sums_tot;

%//cal P_ii_jj
P_ii_jj = squeeze(pair_counts_symm(ii,jj,:))' ./ pair_sums_tot; 

P = P_ii_jj ./ (P_ii.*P_jj);
score1 = -log(P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P(r|ai,aj)P(r)/P(r|ai)/P(r|aj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_ij = sum(pair_counts_symm(ii,jj,:));
P_ii_jj = squeeze(pair_counts_symm(ii,jj,:)) ./ sum_ij;
sum_r_aa = sum(pair_sums_aan,2); % 20
P_ii = pair_sums_aan(ii, :) / sum_r_aa(ii);
P_jj = pair_sums_aan(jj, :) / sum_r_aa(jj);
P_r = pair_sums_tot / totalsum;
P = P_ii_jj' .* P_r ./ (P_ii .* P_jj);
score2 = -log(P);

%shift the end
score1 = score1-score1(cutoff*2+1);
score2 = score2-score2(cutoff*2+1);

%bump
bump_tol = 1;
Nsum_ij = 0;
for i = 1:30
	Nsum_ij = Nsum_ij+pair_counts(ii,jj,i);
	if Nsum_ij>bump_tol
		bump_cut = i;
		break;
	end
end

%bump
if score1(bump_cut+1)>0
	score1(1:bump_cut) = score1(bump_cut+1);
else
	score1(1:bump_cut) = 0;
end
if score2(bump_cut+1)>0
	score2(1:bump_cut) = score2(bump_cut+1);
else
	score2(1:bump_cut) = 0;
end

%x=[2:.1:12];
x=[0.25:0.5:cutoff];
%x=[0.25:0.5:8];
xf=[0.25:0.5:15];
y1 = interp1(xf, score1, x, 'linear'); % // this is the exact pair score rosetta is using
y2 = interp1(xf, score2, x, 'linear');

