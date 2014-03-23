
[list, pair, env, numlist] = load_pair_env_table();

reslist = [
 "ALA"; %1
 "CYS"; %2
 "ASP"; %3
 "GLU"; %4
 "PHE"; %5
 "GLY"; %6
 "HIS"; %7
 "ILE"; %8
 "LYS"; %9
 "LEU"; %10
 "MET"; %11
 "ASN"; %12
 "PRO"; %13
 "GLN"; %14
 "ARG"; %15
 "SER"; %16
 "THR"; %17
 "VAL"; %18
 "TRP"; %19
 "TYR" ];

res_index = zeros(40,1);
env_index = zeros(40,1);
num_res_env = zeros(20,2);
for i = 1:40
  %get aa
  for n = 1:20
    if (list(i,1:3)==reslist(n,:))
      res_index(i)=n;
      break;
    end
  end
  %get env
  env_index(i) = str2num(list(i,5));
  num_res_env(res_index(i),env_index(i))=numlist(i);
end

%split
surf_surf=zeros(20,20,30);
surf_core=zeros(20,20,30);
core_surf=zeros(20,20,30);
core_core=zeros(20,20,30);
for i = 1:40
  ai = res_index(i);
  for j = 1:40
    aj = res_index(j);
    if (env_index(i)==1)
      if (env_index(j)==1)
        %s-s
        surf_surf(ai, aj, :) += pair(i, j, :);
      else
        %s-c
        surf_core(ai, aj, :) += pair(i, j, :);
      end
    else
      if (env_index(j)==1)
        %c-s
        core_surf(ai, aj, :) += pair(i, j, :);
      else
        %c-c
        core_core(ai, aj, :) += pair(i, j, :);
      end
    end
  end
end

%smooth data
L=6;
CONV1 = gausswin(2*L+1, 6);
CONV1 = CONV1/sum(CONV1);
for i = 1:20
  for j= i:20
    smoothed = conv(squeeze(surf_surf(i,j,:)), CONV1);
    surf_surf(i,j,:) = smoothed(L+1:end-L);
    smoothed = conv(squeeze(surf_core(i,j,:)), CONV1);
    surf_core(i,j,:) = smoothed(L+1:end-L);
    smoothed = conv(squeeze(core_surf(i,j,:)), CONV1);
    core_surf(i,j,:) = smoothed(L+1:end-L);
    smoothed = conv(squeeze(core_core(i,j,:)), CONV1);
    core_core(i,j,:) = smoothed(L+1:end-L);
    if (i!=j)
      surf_surf(j,i,:)=surf_surf(i,j,:);
      surf_core(j,i,:)=surf_core(i,j,:);
      core_surf(j,i,:)=core_surf(i,j,:);
      core_core(j,i,:)=core_core(i,j,:);
    end
  end
end

%sum
aa_surf=zeros(20,30);
aa_core=zeros(20,30);
for i = 1:20
  aa_surf(:,:) += reshape((surf_surf(i,:,:)+surf_core(i,:,:)),20,30);
  aa_core(:,:) += reshape((core_surf(i,:,:)+core_core(i,:,:)),20,30);
end

aa_aa  = (surf_surf+surf_core+core_surf+core_core); %20aa to 20aa
pair_counts_symm = aa_aa;
%aa_all = aa_surf+aa_core; %20aa to all


%init
pair_sums_aan = squeeze(sum(pair_counts_symm)); %N(ai, r) = sum_aj N(ai, aj, r)
pair_sums_tot = squeeze(sum(pair_sums_aan)); %N(r) = sum_ai N(ai, r)
totalsum = sum(pair_sums_tot); %N = sum_r N(r)
P_r = pair_sums_tot / totalsum; %P(r) = N(r) / N
sum_r_aa = sum(pair_sums_aan, 2); %N(aa) = sum_r N(ai, r)

%core/surf
sum_r_aa_surf = sum(aa_surf, 2);
sum_r_aa_core = sum(aa_core, 2);

%cal i/j
%if exist("ii", "var")
%  ii
%else
%  ii=10
%end

ii=4

for jj = 1:20
  sum_ij = sum(pair_counts_symm(ii,jj,:));
  P_ii_jj = squeeze(pair_counts_symm(ii,jj,:)) ./ sum_ij;
  %P_ii_surf_func = calc_Pi(aa_surf, aa_core, ii, 1.0);
  %P_ii_core_func = calc_Pi(aa_surf, aa_core, ii, 0.0);
  P_ii = pair_sums_aan(ii, :) / sum_r_aa(ii);
  P_ii_w = calc_Pi(aa_surf, aa_core, ii, num_res_env(jj,1)/sum(num_res_env(jj,:)));
  P_jj = pair_sums_aan(jj, :) / sum_r_aa(jj);
  P_jj_w = calc_Pi(aa_surf, aa_core, jj, num_res_env(ii,1)/sum(num_res_env(ii,:)));
  P = P_ii_jj' .* P_r ./ (P_ii .* P_jj);
  P_w = P_ii_jj' .* P_r ./ (P_ii_w .* P_jj_w);
  score = -log(P);
  score(1:8)=score(9);
  score_w = -log(P_w);
  score_w(1:8)=score_w(9);

  subplot(4,5,jj);
  h=plot(1:30, score, 'r', 1:30, score_w);
  set (h(2), "linewidth", 2)
  axis([0 25 -1 0.6])
  t = reslist(jj,:);
  text(2, score(10), t);
end

