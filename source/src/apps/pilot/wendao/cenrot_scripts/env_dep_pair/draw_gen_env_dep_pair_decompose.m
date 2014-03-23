
[list, pair, env, numlist] = load_pair_env_table();

cutoff=2;
range=12;

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

%mark each table column/row
res_index = zeros(40,1); %residue index
env_index = zeros(40,1); %env: 1-surf, 2-core
num_res_env = zeros(20,2); %residue number for surf & core
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

%split 20x20 x4
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

%add prior dfire-like number
beta = 0.2;
surf_surf += beta*pre_fill_dfire(surf_surf);
surf_core += beta*pre_fill_dfire(surf_core);
core_surf += beta*pre_fill_dfire(core_surf);
core_core += beta*pre_fill_dfire(core_core);

%smooth data
L=5;
CONV1 = gausswin(2*L+1, 4);
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
aa_aa  = (surf_surf+surf_core+core_surf+core_core); %20aa to 20aa
pair_counts_symm = aa_aa;

%init
pair_sums_aan = squeeze(sum(pair_counts_symm)); %N(ai, r) = sum_aj N(ai, aj, r)
pair_sums_tot = squeeze(sum(pair_sums_aan)); %N(r) = sum_ai N(ai, r)
totalsum = sum(pair_sums_tot); %N = sum_r N(r)
P_r = pair_sums_tot / totalsum; %P(r) = N(r) / N
sum_r_aa = sum(pair_sums_aan, 2); %N(aa) = sum_r N(ai, r)

%core/surf
sr_surf_surf = squeeze(sum(surf_surf, 2));
sr_surf_core = squeeze(sum(surf_core, 2));
sr_core_surf = squeeze(sum(core_surf, 2));
sr_core_core = squeeze(sum(core_core, 2));

Sum_r_surf_surf = squeeze(sum(sr_surf_surf,2));
Sum_r_surf_core = squeeze(sum(sr_surf_core,2));
Sum_r_core_surf = squeeze(sum(sr_core_surf,2));
Sum_r_core_core = squeeze(sum(sr_core_core,2));

Nr_surf_surf = sum(sr_surf_surf);
Nr_core_surf = sum(sr_core_surf);
Nr_surf_core = sum(sr_surf_core);
Nr_core_core = sum(sr_core_core);

%pre normalize
%need it or not??? -- yes!
for i = 1:20
  sr_surf_surf(i,:) /= Sum_r_surf_surf(i);
  sr_surf_core(i,:) /= Sum_r_surf_core(i);
  sr_core_surf(i,:) /= Sum_r_core_surf(i);
  sr_core_core(i,:) /= Sum_r_core_core(i);
end
Nr_surf_surf /= sum(Nr_surf_surf);
Nr_core_surf /= sum(Nr_core_surf);
Nr_surf_core /= sum(Nr_surf_core);
Nr_core_core /= sum(Nr_core_core);

%cal i/j
if exist("ii", "var")
  ii
else
  ii=10
end

%fp_pair = fopen("cen_rot_pair.log", "w");

for jj = 4:4
  sum_ij = sum(pair_counts_symm(ii,jj,:));
  P_ii_jj = squeeze(pair_counts_symm(ii,jj,:)) ./ sum_ij;
  P_ii = pair_sums_aan(ii, :) / sum_r_aa(ii);
  P_jj = pair_sums_aan(jj, :) / sum_r_aa(jj);
  %be careful!!
  wi = num_res_env(ii,1)/sum(num_res_env(ii,:));
  wj = num_res_env(jj,1)/sum(num_res_env(jj,:));
  P_ii_w = calc_Pi_sc(sr_surf_surf(ii,:), sr_surf_core(ii,:), sr_core_surf(ii,:), sr_core_core(ii,:), wi, wj);
  P_jj_w = calc_Pi_sc(sr_surf_surf(jj,:), sr_surf_core(jj,:), sr_core_surf(jj,:), sr_core_core(jj,:), wj, wi);
  P = P_ii_jj' .* P_r ./ (P_ii .* P_jj);
  Pr_w = calc_Pi_sc(Nr_surf_surf, Nr_surf_core, Nr_core_surf, Nr_core_core, num_res_env(ii,1)/sum(num_res_env(ii,:)), num_res_env(jj,1)/sum(num_res_env(jj,:)));
  P_w = P_ii_jj' .* P_r ./ (P_ii_w .* P_jj_w);
  P_w2 = P_ii_jj' .* Pr_w ./ (P_ii_w .* P_jj_w);
  P_w3 = P_ii_jj' ./ (P_ii_w .* P_jj_w) .* Pr_w .* Pr_w ./ P_r;

  %score
  score = -log(P);
  score(1:cutoff)=0;
  score_w = -log(P_w);
  score_w(1:cutoff)=0;
  score_w2 = -log(P_w2);
  score_w2(1:cutoff)=0;
  score_w3 = -log(P_w3);
  score_w3(1:cutoff)=0;

  %shift
  %for kk = 1:4
  %  score_w2(kk+20) = score_w2(20)*(4-kk)/4;
  %  score_w(kk+20) = score_w(20)*(4-kk)/4;
  %end
  %score_w2(25:30) = 0;
  %score_w(25:30) = 0;

  %%%%%%%%
  %draw
  %%%%%%%%
  %subplot(4,5,jj);
  %h=plot(1:30, score, 'r', 1:30, score_w, 'b', 1:30, score_w2, 'g', 1:30, score_w3, 'y');
  %h=plot(1:30, score, 'r', 1:30, score_w, 'b', 1:30, score_w2, 'g');
  h=plot(1:30, score, 'r', 1:30, score_w, 'b', 1:30);
  axis([0 25 -1 0.6])
  set (h(1), "linewidth", 2)
  set (h(2), "linewidth", 2)
  %set (h(3), "linewidth", 2)
  t = reslist(jj,:);
  %text(2, 0, t);

  %%%%%%%%%%%%
  %% output
  %%%%%%%%%%%%
  %y2 = score_w2(1:range*2);
  %fprintf(fp_pair, "PAIR: %s %s ", reslist(ii,:), reslist(jj,:));
  %fprintf(fp_pair, "%.2f %.2f ", 0.25, 0.25+(range*2-1)*0.5);
  %fprintf(fp_pair, "%d ", length(y2));
  %for k=1:length(y2)
  %  fprintf(fp_pair, "%f ", y2(k));
  %end
  %fprintf(fp_pair, "\n");
end
%fprintf(fp_pair, "# END\n");
%fclose(fp_pair);

%% debug
%ii=4
%jj=8
%subplot(2,2,1)
%plot(1:30,squeeze(sr_surf_surf(ii,:)), 1:30, squeeze(sr_surf_surf(jj,:)))
%subplot(2,2,2)
%plot(1:30,squeeze(sr_surf_core(ii,:)), 1:30, squeeze(sr_core_surf(jj,:)))
%subplot(2,2,3)
%plot(1:30,squeeze(sr_core_surf(ii,:)), 1:30, squeeze(sr_surf_core(jj,:)))
%subplot(2,2,4)
%plot(1:30,squeeze(sr_core_core(ii,:)), 1:30, squeeze(sr_core_core(jj,:)))

