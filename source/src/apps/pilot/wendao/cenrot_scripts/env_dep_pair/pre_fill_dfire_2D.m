function [extra] = pre_fill_dfire_2D(raw, NAA)
%alpha = 1.61;
alpha = 2.0;

Rbin = 0.5;
Rcut = 14.75;
r = linspace(Rbin/2, Rcut, 30);
sumr = squeeze(sum(raw,4));
Ncut = squeeze(sumr(:,:,30));
extra = zeros(NAA, NAA, 30, 16);

for i = 1:NAA
  for j = 1:NAA
    for k = 1:30
      extra( i, j, k, : ) = (r(k)/Rcut)^alpha * Ncut(i, j) / 16;
    end
  end
end

