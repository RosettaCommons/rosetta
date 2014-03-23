
[list, th1_table, th2_table, dih_table, numlist] = load_angle_r_table();

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

%buildlist
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
surf_surf_th1=zeros(20,20,30,16);
surf_core_th1=zeros(20,20,30,16);
core_surf_th1=zeros(20,20,30,16);
core_core_th1=zeros(20,20,30,16);
surf_surf_th2=zeros(20,20,30,16);
surf_core_th2=zeros(20,20,30,16);
core_surf_th2=zeros(20,20,30,16);
core_core_th2=zeros(20,20,30,16);
for i = 1:40
  ai = res_index(i);
  for j = 1:40
    aj = res_index(j);
    if (env_index(i)==1)
      if (env_index(j)==1)
        %s-s
        surf_surf_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        surf_surf_th2(ai, aj, :, :) += th2_table(i, j, :, :);
      else
        %s-c
        surf_core_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        surf_core_th2(ai, aj, :, :) += th2_table(i, j, :, :);
      end
    else
      if (env_index(j)==1)
        %c-s
        core_surf_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        core_surf_th2(ai, aj, :, :) += th2_table(i, j, :, :);
      else
        %c-c
        core_core_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        core_core_th2(ai, aj, :, :) += th2_table(i, j, :, :);
      end
    end
  end
end

%smooth

%sum
aa_aa_th1 = (surf_surf_th1 + surf_core_th1 + core_surf_th1 + core_core_th1)+0.1; %20aa to 20aa
aa_aa_th2 = (surf_surf_th2 + surf_core_th2 + core_surf_th2 + core_core_th2)+0.1; %20aa to 20aa

aa_surf_da1=zeros(20,30,16); %aa theta1 -> surf
aa_core_da1=zeros(20,30,16); %aa theta1 -> core
aa_surf_da2=zeros(20,30,16); %surf theta1 -> aa
aa_core_da2=zeros(20,30,16); %core theta1 -> aa
for i = 1:20
  for j = 1:20
    aa_surf_da1(i,:,:) += reshape(squeeze(surf_surf_th1(i,j,:,:)+core_surf_th1(i,j,:,:)),1,30,16);
    aa_surf_da1(i,:,:) += reshape(squeeze(surf_core_th2(j,i,:,:)+surf_surf_th2(j,i,:,:)),1,30,16);
    aa_core_da1(i,:,:) += reshape(squeeze(surf_core_th1(i,j,:,:)+core_core_th1(i,j,:,:)),1,30,16);
    aa_core_da1(i,:,:) += reshape(squeeze(core_core_th2(j,i,:,:)+core_surf_th2(j,i,:,:)),1,30,16);
    aa_surf_da2(i,:,:) += reshape(squeeze(surf_surf_th1(j,i,:,:)+surf_core_th1(j,i,:,:)),1,30,16);
    aa_surf_da2(i,:,:) += reshape(squeeze(core_surf_th2(i,j,:,:)+surf_surf_th2(i,j,:,:)),1,30,16);
    aa_core_da2(i,:,:) += reshape(squeeze(core_surf_th1(j,i,:,:)+core_core_th1(j,i,:,:)),1,30,16);
    aa_core_da2(i,:,:) += reshape(squeeze(core_core_th2(i,j,:,:)+surf_core_th2(i,j,:,:)),1,30,16);
  end
end

%cal
ntype = 20;
dat_ang = aa_aa_th1;
dat_a_a_r_ang = zeros(ntype, ntype, 30, 16);
for aa = 1:ntype
  for bb = 1:ntype
    dat_ang(aa,bb,:,:) += aa_aa_th2(bb,aa,:,:);
    for kk = 1:30
      dat_ang(aa,bb,kk,:) = squeeze(dat_ang(aa,bb,kk,:));
    end
    %don't smooth with zero
    dat_a_a_r_ang(aa,bb,:,:) = dat_ang(aa,bb,:,:);
    %dat_a_a_r_ang(aa,bb,cc:30,:) = smooth_2D(dat_ang(aa,bb,cc:30,:), [15, 0.75, 2], [180, 9, 2]);
  end
end

sum_a_a_r = sum(dat_a_a_r_ang,4); %sum over ang
dat_a_r_ang1 = squeeze(sum(dat_a_a_r_ang, 2)); %sum over the second aa
sum_a1_r = sum(dat_a_r_ang1, 3); %sum over ang
dat_a_r_ang2 = squeeze(sum(dat_a_a_r_ang)); %sum over the first aa
sum_a2_r = sum(dat_a_r_ang2, 3); %sum over ang

dat_r_ang = squeeze(sum(dat_a_r_ang1)); %sum over all aa
sum_r = sum(dat_r_ang,2); %sum over ang

p_a_a_r_ang = zeros(ntype, ntype, 30, 16);
p_a_r_ang1 = zeros(ntype, 30, 16);
p_a_r_ang2 = zeros(ntype, 30, 16);
p_r_ang = zeros(30, 16);

aa_surf_dat_r_ang = squeeze(sum(aa_surf));
aa_core_dat_r_ang = squeeze(sum(aa_core));
aa_surf_sum_r = sum(aa_surf_dat_r_ang, 2);
aa_core_sum_r = sum(aa_core_dat_r_ang, 2);
aa_surf_p_r_ang = zeros(30,16);
aa_core_p_r_ang = zeros(30,16);

for kk = 1:30
  for ii = 1:ntype
    for jj = 1:ntype
      p_a_a_r_ang(ii,jj,kk,:) = dat_a_a_r_ang(ii,jj,kk,:) / sum_a_a_r(ii,jj,kk);
    end
    p_a_r_ang1(ii,kk,:) = dat_a_r_ang1(ii,kk,:) / sum_a1_r(ii,kk);
    p_a_r_ang2(ii,kk,:) = dat_a_r_ang2(ii,kk,:) / sum_a2_r(ii,kk);
  end
  p_r_ang(kk,:) = dat_r_ang(kk,:) / sum_r(kk);
  aa_surf_p_r_ang(kk,:) = aa_surf_dat_r_ang(kk,:) / aa_surf_sum_r(kk);
  aa_core_p_r_ang(kk,:) = aa_core_dat_r_ang(kk,:) / aa_core_sum_r(kk);
end

if exist("iii", "var")
  iii
else
  iii=4
end
if exist("jjj", "var")
  jjj
else
  jjj=8
end

p_ii_jj = squeeze(p_a_a_r_ang( iii, jjj, :, : ));
p_ii = p_a_r_ang1(iii, :, :);
%p_ii_w = calc_Pi_ang(aa_surf_da1, aa_core_da1, iii, num_res_env(jjj,1)/sum(num_res_env(jjj,:)));
p_ii_w = calc
p_jj = p_a_r_ang2(jj, :, :);
%p_jj_w = calc_Pi_ang(aa_surf_da2, aa_core_da2, jjj, num_res_env(iii,1)/sum(num_res_env(iii,:)));
p_ii_p_jj = squeeze(p_ii .* p_jj);
p_ii_p_jj_w = squeeze(p_ii_w .* p_jj_w);
p = (p_ii_jj .* p_r_ang) ./ (p_ii_p_jj);
p_w = (p_ii_jj .* p_r_ang) ./ (p_ii_p_jj_w);
score = -log(p);
score_w = -log(p_w);
score(score>1) = 1.0;
score(score<-1) = -1.0;
score_w(score_w>1) = 1.0;
score_w(score_w<-1) = -1.0;

%printf("score(%d,%d): max=%f, min=%f\n", ii, jj, max(max(score)), min(min(score)));

[xx, yy] = meshgrid(1:16, 1:30);
%contourf(xx,yy,score);
%plot(1:16, p_ii_w(30,:), 1:16, squeeze(p_ii(1,30,:)),1:16,squeeze(aa_surf_da1(4,30,:))/sum(squeeze(aa_surf_da1(4,30,:))),1:16,squeeze(aa_core_da1(4,30,:))/sum(squeeze(aa_core_da1(4,30,:))))
plot(1:16,score(30,:),"r",1:16,score_w(30,:),"b")
