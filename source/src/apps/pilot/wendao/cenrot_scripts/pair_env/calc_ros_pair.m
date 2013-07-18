function y=calc_ros_pair(pair,ii,jj)

x=[3:.1:12];
x2 = x.^2;

% // r++
cen_dist5_pad = 0.5; cen_dist6_pad = 0.6; cen_dist7_pad = 0.65; cen_dist10_pad = 1.0; cen_dist12_pad = 1.2;
cen_dist5_pad_plus = cen_dist5_pad  + 25.0;
cen_dist6_pad_plus = cen_dist6_pad + 36.0;
cen_dist7_pad_plus = cen_dist7_pad  + 56.25;
cen_dist10_pad_plus = cen_dist10_pad + 100.0;
cen_dist12_pad_plus = cen_dist12_pad + 144.0;
cen_dist5_pad_minus = cen_dist5_pad  - 25.0;
cen_dist7_pad_minus = cen_dist7_pad  - 56.25;
cen_dist10_pad_minus = cen_dist10_pad - 100.0;
cen_dist12_pad_minus = cen_dist12_pad - 144.0;
cen_dist5_pad_hinv = 0.5 / cen_dist5_pad;
cen_dist6_pad_hinv = 0.5 / cen_dist6_pad;
cen_dist7_pad_hinv =  0.5 / cen_dist7_pad;
cen_dist10_pad_hinv = 0.5 / cen_dist10_pad;
cen_dist12_pad_hinv = 0.5 / cen_dist12_pad;

x0=x2;
x0( x2>=cen_dist10_pad_plus ) = 4+max( ( x2( x2>=cen_dist10_pad_plus )+cen_dist12_pad_minus ) * cen_dist12_pad_hinv, 0);
x0( x2>=cen_dist7_pad_plus & x2< cen_dist10_pad_plus) = ...
	3 + max( ( x2( x2>=cen_dist7_pad_plus & x2< cen_dist10_pad_plus) + cen_dist10_pad_minus ) * cen_dist10_pad_hinv, 0);
x0( x2>=cen_dist5_pad_plus & x2< cen_dist7_pad_plus) = ...
	2 + max( ( x2( x2>=cen_dist5_pad_plus & x2< cen_dist7_pad_plus)  + cen_dist7_pad_minus ) * cen_dist7_pad_hinv , 0);
x0( x2< cen_dist5_pad_plus) = ...
	1 + max( ( x2( x2< cen_dist5_pad_plus) + cen_dist5_pad_minus ) * cen_dist5_pad_hinv, 0);

y = interp1(1:5,pair(:,ii,jj),x0,'linear'); % // this is the exact pair score rosetta is using
