
[list, pair, env, numlist] = load_pair_env_table();

reslist = [
 "LYS"; %1
 "GLU"; %2
 "ASP"; %3
 "ARG"; %4
 "GLN"; %5
 "ASN"; %6
 "PRO"; %7
 "SER"; %8
 "THR"; %9
 "HIS"; %10
 "GLY"; %11
 "ALA"; %12
 "TYR"; %13
 "MET"; %14
 "TRP"; %15
 "VAL"; %16
 "LEU"; %17
 "PHE"; %18
 "ILE"; %19
 "CYS" ];

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

for i = 1:20
  printf("%s: %f\n", reslist(i,:), num_res_env(i,1)/(sum(num_res_env(i,:))))
end

surf_surf_close = sum(surf_surf(:,:,1:14),3);
surf_core_close = sum(surf_core(:,:,1:14),3);
core_surf_close = sum(core_surf(:,:,1:14),3);
core_core_close = sum(core_core(:,:,1:14),3);

[xx, yy] = meshgrid(1:20, 1:20);
mat = calc_pair_matrix(core_core_close);
subplot(2,2,1);
imshow(1-mat/2)
mat = calc_pair_matrix(surf_surf_close);
subplot(2,2,2);
imshow(1-mat/2)
mat = calc_pair_matrix(surf_core_close+core_surf_close);
subplot(2,2,3);
imshow(1-mat/2)
mat = calc_pair_matrix(surf_surf_close+surf_core_close+core_surf_close);
subplot(2,2,4);
imshow(1-mat/2)

