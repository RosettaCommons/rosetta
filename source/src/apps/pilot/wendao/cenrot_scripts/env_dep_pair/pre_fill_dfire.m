function [extra] = pre_fill_dfire(raw)
alpha = 1.61;
%alpha = 2.0;

Rbin = 0.5;
Rcut = 14.75;
r = linspace(Rbin/2, Rcut, 30);
Ncut = squeeze(raw(:,:,30));
extra = zeros(20, 20, 30);

for i = 1:20
  for j = 1:20
    for k = 1:30
      extra( i, j, k ) = (r(k)/Rcut)^alpha * Ncut(i, j);
    end
  end
end
