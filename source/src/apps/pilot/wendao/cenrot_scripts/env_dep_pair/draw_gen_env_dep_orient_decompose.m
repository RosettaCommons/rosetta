
[list, th1_table, th2_table, dih_table, numlist] = load_angle_r_table();

debug = 1;

reslist = [
 "TYR"; %1
 "CYS"; %2
 "ASP"; %3
 "GLU"; %4
 "PHE"; %5
 "TRP"; %6
 "HIS"; %7
 "ILE"; %8
 "LYS"; %9
 "LEU"; %10
 "MET"; %11
 "ASN"; %12
 "VAL"; %13
 "GLN"; %14
 "ARG"; %15
 "SER"; %16
 "THR"; %17
 "PRO"; %18
 "ALA"; %19
 "GLY" ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% should not contain ALA GLY PRO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAA = 17;
Rcut = 12;
Bcut = 24;

%buildlist
res_index = zeros(40,1);
env_index = zeros(40,1);
num_res_env = zeros(NAA,2);
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
surf_surf_th1=zeros(NAA,NAA,30,16);
surf_core_th1=zeros(NAA,NAA,30,16);
core_surf_th1=zeros(NAA,NAA,30,16);
core_core_th1=zeros(NAA,NAA,30,16);
for i = 1:40
  ai = res_index(i);
  if (ai>NAA)
    continue
  end
  for j = 1:40
    aj = res_index(j);
    if (aj>NAA)
      continue
    end
    if (env_index(i)==1)
      if (env_index(j)==1)
        %s-s
        surf_surf_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        surf_surf_th1(ai, aj, :, :) += th2_table(j, i, :, :);
      else
        %s-c
        surf_core_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        surf_core_th1(ai, aj, :, :) += th2_table(j, i, :, :);
      end
    else
      if (env_index(j)==1)
        %c-s
        core_surf_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        core_surf_th1(ai, aj, :, :) += th2_table(j, i, :, :);
      else
        %c-c
        core_core_th1(ai, aj, :, :) += th1_table(i, j, :, :);
        core_core_th1(ai, aj, :, :) += th2_table(j, i, :, :);
      end
    end
  end
end

%add prior-data
beta = 0.2;
surf_surf_th1 += beta*pre_fill_dfire_2D(surf_surf_th1, NAA);
surf_core_th1 += beta*pre_fill_dfire_2D(surf_core_th1, NAA);
core_surf_th1 += beta*pre_fill_dfire_2D(core_surf_th1, NAA);
core_core_th1 += beta*pre_fill_dfire_2D(core_core_th1, NAA);

%smooth data
%for i = 1:NAA
%  for j = i:NAA
%    smoothed = smooth_2D(surf_surf_th1(i,j,:,:), [15, 0.1, 2], [360, 40, 2]);
%    surf_surf_th1(i,j,:,:) = smoothed;
%    smoothed = smooth_2D(surf_core_th1(i,j,:,:), [15, 0.1, 2], [360, 40, 2]);
%    surf_core_th1(i,j,:,:) = smoothed;
%    smoothed = smooth_2D(core_surf_th1(i,j,:,:), [15, 0.1, 2], [360, 40, 2]);
%    core_surf_th1(i,j,:,:) = smoothed;
%    smoothed = smooth_2D(core_core_th1(i,j,:,:), [15, 0.1, 2], [360, 40, 2]);
%    core_core_th1(i,j,:,:) = smoothed;
%    if (i!=j)
%      surf_surf_th1(j,i,:,:) = surf_surf_th1(i,j,:,:);
%      surf_core_th1(j,i,:,:) = surf_core_th1(i,j,:,:);
%      core_surf_th1(j,i,:,:) = core_surf_th1(i,j,:,:);
%      core_core_th1(j,i,:,:) = core_core_th1(i,j,:,:);
%    end
%  end
%end

%sum symm
aasurf_surf_da1=zeros(NAA,30,16); %surf aa theta1 -> surf
aacore_surf_da1=zeros(NAA,30,16); %core aa theta1 -> surf
aasurf_core_da1=zeros(NAA,30,16); %surf aa theta1 -> core
aacore_core_da1=zeros(NAA,30,16); %core aa theta1 -> core
aasurf_surf_da2=zeros(NAA,30,16); %surf aa -> surf theta2
aacore_surf_da2=zeros(NAA,30,16); %core aa -> surf theta2
aasurf_core_da2=zeros(NAA,30,16); %surf aa -> core theta2
aacore_core_da2=zeros(NAA,30,16); %core aa -> core theta2

aasurf_surf_da1=squeeze(sum(surf_surf_th1,2));
aacore_surf_da1=squeeze(sum(core_surf_th1,2));
aasurf_core_da1=squeeze(sum(surf_core_th1,2));
aacore_core_da1=squeeze(sum(core_core_th1,2));
aasurf_surf_da2=squeeze(sum(surf_surf_th1));
aacore_surf_da2=squeeze(sum(surf_core_th1));
aasurf_core_da2=squeeze(sum(core_surf_th1));
aacore_core_da2=squeeze(sum(core_core_th1));

Na1_r_surf_surf = squeeze(sum(aasurf_surf_da1));
Na1_r_surf_core = squeeze(sum(aasurf_core_da1));
Na1_r_core_surf = squeeze(sum(aacore_surf_da1));
Na1_r_core_core = squeeze(sum(aacore_core_da1));
Na2_r_surf_surf = squeeze(sum(aasurf_surf_da2));
Na2_r_surf_core = squeeze(sum(aasurf_core_da2));
Na2_r_core_surf = squeeze(sum(aacore_surf_da2));
Na2_r_core_core = squeeze(sum(aacore_core_da2));

%%pre-normalize
for i = 1:NAA
  aasurf_surf_da1(i,:,:) = calc_NtoP(aasurf_surf_da1(i,:,:));
  aasurf_core_da1(i,:,:) = calc_NtoP(aasurf_core_da1(i,:,:));
  aacore_surf_da1(i,:,:) = calc_NtoP(aacore_surf_da1(i,:,:));
  aacore_core_da1(i,:,:) = calc_NtoP(aacore_core_da1(i,:,:));
  aasurf_surf_da2(i,:,:) = calc_NtoP(aasurf_surf_da2(i,:,:));
  aacore_surf_da2(i,:,:) = calc_NtoP(aacore_surf_da2(i,:,:));
  aasurf_core_da2(i,:,:) = calc_NtoP(aasurf_core_da2(i,:,:));
  aacore_core_da2(i,:,:) = calc_NtoP(aacore_core_da2(i,:,:));
end
Na1_r_surf_surf = calc_NtoP(Na1_r_surf_surf);
Na1_r_surf_core = calc_NtoP(Na1_r_surf_core);
Na1_r_core_surf = calc_NtoP(Na1_r_core_surf);
Na1_r_core_core = calc_NtoP(Na1_r_core_core);
Na2_r_surf_surf = calc_NtoP(Na2_r_surf_surf);
Na2_r_surf_core = calc_NtoP(Na2_r_surf_core);
Na2_r_core_surf = calc_NtoP(Na2_r_core_surf);
Na2_r_core_core = calc_NtoP(Na2_r_core_core);

aa_aa_th1 = (surf_surf_th1 + surf_core_th1 + core_surf_th1 + core_core_th1); %NAAaa to NAAaa

%cal
ntype = NAA;
dat_a_a_r_ang = aa_aa_th1;

%old way
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

for kk = 1:30
  for ii = 1:ntype
    for jj = 1:ntype
      p_a_a_r_ang(ii,jj,kk,:) = dat_a_a_r_ang(ii,jj,kk,:) / sum_a_a_r(ii,jj,kk);
    end
    p_a_r_ang1(ii,kk,:) = dat_a_r_ang1(ii,kk,:) / sum_a1_r(ii,kk);
    p_a_r_ang2(ii,kk,:) = dat_a_r_ang2(ii,kk,:) / sum_a2_r(ii,kk);
  end
  p_r_ang(kk,:) = dat_r_ang(kk,:) / sum_r(kk);
end

%shift
shift = zeros(30,16);
%for i = 20:30
%  shift(i,:) = (i-20)/10*0.04;
%end

%debug
if debug
  f = figure();
  set(f, "visible", "off")
else
  fp_ang = fopen("cen_rot_pair_ang.log", "w");
end

%if exist("iii", "var")
%  iii
%else
%  iii=4
%end
%if exist("jjj", "var")
%  jjj
%else
%  jjj=8
%end

%for iii = 1:NAA
%for jjj = 1:NAA

iii = 4
jjj = 13

wi = num_res_env(iii,1)/sum(num_res_env(iii,:));
wj = num_res_env(jjj,1)/sum(num_res_env(jjj,:));
p_ii_jj = squeeze(p_a_a_r_ang( iii, jjj, :, : ));
p_ii = p_a_r_ang1(iii, :, :);
p_ii_w = calc_Pi_ang_sc(aasurf_surf_da1(iii,:,:), aasurf_core_da1(iii,:,:), aacore_surf_da1(iii,:,:), aacore_core_da1(iii,:,:), wi, wj);
p_jj = p_a_r_ang2(jjj, :, :);
p_jj_w = calc_Pi_ang_sc(aasurf_surf_da2(jjj,:,:), aasurf_core_da2(jjj,:,:), aacore_surf_da2(jjj,:,:), aacore_core_da2(jjj,:,:), wj, wi);
p_ii_p_jj = squeeze(p_ii .* p_jj);
p_ii_p_jj_w = squeeze(p_ii_w .* p_jj_w);
p = (p_ii_jj .* p_r_ang) ./ (p_ii_p_jj);
p_r_ang_w = calc_Pi_ang_sc(Na1_r_surf_surf, Na1_r_surf_core, Na1_r_core_surf, Na1_r_core_core, wi, wj);
p_w = (p_ii_jj .* p_r_ang_w + shift) ./ (p_ii_p_jj_w + shift);
p_w0 = (p_ii_jj .* p_r_ang + shift) ./ (p_ii_p_jj_w + shift);
score0 = -log(p);
score_w = -log(p_w);
score_w0 = -log(p_w0);
score_w = smooth_2D(score_w, [15,1.5,2], [360,36,2]);
score0 = smooth_2D(score0, [15,1.5,2], [360,36,2]);

score = score_w(1:Bcut,:);

%debug
if debug
  printf("score_w(%d,%d): max=%f, min=%f\n", iii, jjj, max(max(score)), min(min(score)));
  [xx, yy] = meshgrid(1:16, 1:Bcut);
  contourf(xx,yy,score);
  fname = [ "score_", int2str(iii), "_", int2str(jjj), ".png" ];
  print(fname, "-dpng")
else
  fprintf(fp_ang, "ANGLE1: %s %s -0.25 12.25 26 -1.0625 1.0625 18\n", reslist(iii,:), reslist(jjj,:));
  for yy=1:18
    fprintf(fp_ang, "%.6f ", 0);
  end
  fprintf(fp_ang, "\n");
  for xx=1:Bcut
    fprintf(fp_ang, "%.6f ", score(xx,1));
    for yy=1:16
      fprintf(fp_ang, "%.6f ", score(xx,yy));
    end
    fprintf(fp_ang, "%.6f ", score(xx,16));
    fprintf(fp_ang, "\n");
  end
  for yy=1:18
    fprintf(fp_ang, "%.6f ", 0);
  end
  fprintf(fp_ang, "\n");
end

%end
%end

if debug
  printf("Done!\n");
  set(f, "visible", "on")
else
  fclose(fp_ang);
end

%% debug %%

%subplot(2,1,1)
%[xx, yy] = meshgrid(1:16, 1:30);
%contourf(xx,yy,score_w);
%plot(1:16, p_ii_w(30,:), 1:16, squeeze(p_ii(1,30,:)),1:16,squeeze(aa_surf_da1(4,30,:))/sum(squeeze(aa_surf_da1(4,30,:))),1:16,squeeze(aa_core_da1(4,30,:))/sum(squeeze(aa_core_da1(4,30,:))))
%subplot(2,1,2)
%plot(1:16,score(24,:),"r",1:16,score_w(24,:),"b")
%axis([1 16 -0.2 0.2])
%plot(1:16,p_ii_jj(30,:),"r",1:16,p_r_ang_w(30,:),"g",1:16, Na1_r_surf_core(30,:)/sum(Na1_r_surf_core(30,:)), 'c', 1:16, p_ii_w(30,:), 'm', 1:16, p_jj_w(30,:), 'b', 1:16, p_r_ang(30,:))

%subplot(2,2,1)
%plot(1:16,squeeze(aasurf_surf_da1(iii,30,:)), 1:16, squeeze(aasurf_surf_da2(jjj,30,:)))
%subplot(2,2,2)
%plot(1:16,squeeze(aasurf_core_da1(iii,30,:)), 1:16, squeeze(aacore_surf_da2(jjj,30,:)))
%subplot(2,2,3)
%plot(1:16,squeeze(aacore_surf_da1(iii,30,:)), 1:16, squeeze(aasurf_core_da2(jjj,30,:)))
%subplot(2,2,4)
%plot(1:16,squeeze(aacore_core_da1(iii,30,:)), 1:16, squeeze(aacore_core_da2(jjj,30,:)))

%subplot(4,2,1);
%contourf(xx,yy, Na1_r_surf_surf )
%subplot(4,2,2);
%contourf(xx,yy, Na2_r_surf_surf )
%subplot(4,2,3);
%contourf(xx,yy, Na1_r_core_surf )
%subplot(4,2,4);
%contourf(xx,yy, Na2_r_surf_core )
%subplot(4,2,5);
%contourf(xx,yy, Na1_r_surf_core )
%subplot(4,2,6);
%contourf(xx,yy, Na2_r_core_surf )
%subplot(4,2,7);
%contourf(xx,yy, Na1_r_core_core )
%subplot(4,2,8);
%contourf(xx,yy, Na2_r_core_core )

