
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
aa_surf=zeros(20,30,16);
aa_core=zeros(20,30,16);
surf_surf=zeros(20,30,16);
surf_core=zeros(20,30,16);
core_surf=zeros(20,30,16);
core_core=zeros(20,30,16);
for i = 1:20
  aa_surf += reshape((surf_surf_th1(i,:,:,:)+surf_core_th1(i,:,:,:)),20,30,16)+0.05;
  aa_surf += reshape((core_surf_th2(:,i,:,:)+surf_surf_th2(:,i,:,:)),20,30,16)+0.05;
  aa_core += reshape((core_surf_th1(i,:,:,:)+core_core_th1(i,:,:,:)),20,30,16)+0.05;
  aa_core += reshape((core_core_th2(:,i,:,:)+surf_core_th2(:,i,:,:)),20,30,16)+0.05;

  surf_surf += squeeze(surf_surf_th1(i,:,:,:))+squeeze(surf_surf_th2(:,i,:,:))+0.05;
  surf_core += squeeze(surf_core_th1(i,:,:,:))+squeeze(core_surf_th2(:,i,:,:))+0.05;
  core_surf += squeeze(core_surf_th1(i,:,:,:))+squeeze(surf_core_th2(:,i,:,:))+0.05;
  core_core += squeeze(core_core_th1(i,:,:,:))+squeeze(core_core_th2(:,i,:,:))+0.05;
end
aa_aa_th1 = (surf_surf_th1 + surf_core_th1 + core_surf_th1 + core_core_th1)+0.1; %20aa to 20aa
aa_aa_th2 = (surf_surf_th2 + surf_core_th2 + core_surf_th2 + core_core_th2)+0.1; %20aa to 20aa

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
surf_surf_dat_r_ang = squeeze(sum(surf_surf));
surf_core_dat_r_ang = squeeze(sum(surf_core));
core_surf_dat_r_ang = squeeze(sum(core_surf));
core_core_dat_r_ang = squeeze(sum(core_core));
aa_surf_sum_r = sum(aa_surf_dat_r_ang, 2);
aa_core_sum_r = sum(aa_core_dat_r_ang, 2);
surf_surf_sum_r = sum(surf_surf_dat_r_ang, 2);
surf_core_sum_r = sum(surf_core_dat_r_ang, 2);
core_surf_sum_r = sum(core_surf_dat_r_ang, 2);
core_core_sum_r = sum(core_core_dat_r_ang, 2);
aa_surf_p_r_ang = zeros(30,16);
aa_core_p_r_ang = zeros(30,16);
surf_surf_p_r_ang = zeros(30,16);
surf_core_p_r_ang = zeros(30,16);
core_surf_p_r_ang = zeros(30,16);
core_core_p_r_ang = zeros(30,16);

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
  surf_surf_p_r_ang(kk,:) = surf_surf_dat_r_ang(kk,:) / surf_surf_sum_r(kk);
  surf_core_p_r_ang(kk,:) = surf_core_dat_r_ang(kk,:) / surf_core_sum_r(kk);
  core_surf_p_r_ang(kk,:) = core_surf_dat_r_ang(kk,:) / core_surf_sum_r(kk);
  core_core_p_r_ang(kk,:) = core_core_dat_r_ang(kk,:) / core_core_sum_r(kk);
end

[xx, yy] = meshgrid(1:16, 1:30);
%contourf(xx,yy,score);

SS = smooth_2D(surf_surf_p_r_ang, [15,1,2], [2,0.1,2]);
SC = smooth_2D(surf_core_p_r_ang, [15,1,2], [2,0.1,2]);
CS = smooth_2D(core_surf_p_r_ang, [15,1,2], [2,0.1,2]);
CC = smooth_2D(core_core_p_r_ang, [15,1,2], [2,0.1,2]);

