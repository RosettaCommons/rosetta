
[list, dat_ang1, dat_ang2, dat_dih] = load_angle_r_table();

debug = 0;
cutoff = 10; %bin tol
dis_cut = 24; %pair term cutoff 12A
range = 0.6;
ntype = size(list)(1);

dat_ang = dat_ang1;
dat_a_a_r_ang = zeros(ntype, ntype, 30, 16);
bump = zeros(ntype, ntype);
for aa = 1:ntype
	for bb = 1:ntype
		dat_ang(aa,bb,:,:) += dat_ang2(bb,aa,:,:);
		for kk = 1:30
			dat_ang(aa,bb,kk,:) = squeeze(dat_ang(aa,bb,kk,:));
		end
		tmpdata = squeeze(dat_ang(aa,bb,:,:));
		tmpsumr = sum(tmpdata, 2);
		for cc = 1:30
			if tmpsumr(cc)>cutoff*16
				bump(aa,bb) = cc;
				break;
			end
		end
		%don't smooth with zero
		dat_a_a_r_ang(aa,bb,cc:30,:) = dat_ang(aa,bb,cc:30,:);
		%dat_a_a_r_ang(aa,bb,cc:30,:) = smooth_2D(dat_ang(aa,bb,cc:30,:), [15, 0.75, 2], [180, 9, 2]);
	end
end

sum_a_a_r = sum(dat_a_a_r_ang,4); %sum over ang
dat_a_r_ang1 = squeeze(sum(dat_a_a_r_ang, 2)); %sum over the second aa
sum_a1_r = sum(dat_a_r_ang1, 3); %sum over ang
dat_a_r_ang2 = squeeze(sum(dat_a_a_r_ang)); %sum over the first aa
sum_a2_r = sum(dat_a_r_ang2, 3); %sum over ang

%dat_a_r_ang1 and 2 should give the save answer
dat_r_ang = squeeze(sum(dat_a_r_ang1)); %sum over all aa
sum_r = sum(dat_r_ang,2); %sum over ang

p_a_a_r_ang = zeros(ntype, ntype, 30, 16);
p_a_r_ang1 = zeros(ntype, 30, 16);
p_a_r_ang2 = zeros(ntype, 30, 16);
p_r_ang = zeros(30, 16);
for kk = 1:30
	for ii = 1:ntype
		for jj = 1:ntype
			if sum_a_a_r(ii,jj,kk)>cutoff*16
				p_a_a_r_ang(ii,jj,kk,:) = dat_a_a_r_ang(ii,jj,kk,:) / sum_a_a_r(ii,jj,kk);
			else
				p_a_a_r_ang(ii,jj,kk,:) = 1.0/16;
			end
		end
		if sum_a1_r(ii,kk)>cutoff*16*ntype
			p_a_r_ang1(ii,kk,:) = dat_a_r_ang1(ii,kk,:) / sum_a1_r(ii,kk);
		else
			p_a_r_ang1(ii,kk,:) = 1.0/16;
		end
		if sum_a2_r(ii,kk)>cutoff*16*ntype
			p_a_r_ang2(ii,kk,:) = dat_a_r_ang2(ii,kk,:) / sum_a2_r(ii,kk);
		else
			p_a_r_ang2(ii,kk,:) = 1.0/16;
		end
	end
	if sum_r(kk)>cutoff*16*ntype*ntype
		p_r_ang(kk,:) = dat_r_ang(kk,:) / sum_r(kk);
	else
		p_r_ang(kk,:) = 1.0/16;
	end
end

p_r_ang = smooth_2D(p_r_ang, [15, 1, 2], [2, 0.2, 2]);

%make up a prior count to force the score down to 0
%prior = zeros(30,16);
%h = 0.04;
%shift = 16;
%for i = shift:30
%	prior(i,:) = -(i-30)*(i-30)/(30-shift)/(30-shift)*h+h;
%end
shift = 15;
prior = zeros(30,16);
for i = shift:30
	prior(i,:) = (i-shift)*0.03/(30-shift);
end

%debug
if debug
	[xx, yy] = meshgrid(1:16, 1:24);
	f = figure();
	set(f, "visible", "off")
else
	fp_ang = fopen("cen_rot_pair_ang.log", "w");
end

for ii = 1:ntype
	for jj = 1:ntype
		%skip GLY/ALA
		if bump(ii,jj)==0
			continue;
		end
		
		mask0 = squeeze(dat_ang(ii,jj,:,:))<cutoff; %no enough data
		p_ii_jj = ones(30, 16)/16;
		p_ii_jj(bump(ii,jj):30,:) = squeeze(p_a_a_r_ang(ii, jj, bump(ii,jj):30, :));
		mask1 = p_ii_jj==0;
		p_ii_jj(bump(ii,jj):30,:) = smooth_2D(p_ii_jj(bump(ii,jj):30,:),[15, 1, 2], [2, 0.2, 2]); 
		p_ii = ones(30, 16)/16;
		p_jj = ones(30, 16)/16;
		p_ii(bump(ii,jj):30,:) = p_a_r_ang1(ii,bump(ii,jj):30,:); 
		p_jj(bump(ii,jj):30,:) = p_a_r_ang2(jj,bump(ii,jj):30,:);
		p_ii_p_jj = p_ii .* p_jj;
		mask2 = p_ii_p_jj==0;
		p_ii_p_jj(bump(ii,jj):30,:) = smooth_2D(p_ii_p_jj(bump(ii,jj):30,:),[15, 1, 2], [2, 0.2, 2]); 
		
		mask = (mask0+mask1+mask2)>0;
		p_ii_jj(mask)=1.0/16;
		p_ii_p_jj(mask)=1.0/16/16;
		
		p = (p_ii_jj .* p_r_ang + prior) ./ (p_ii_p_jj + prior);
		p(mask) = 1;
		p(1:bump(ii,jj)-1,:) = 1;
		score = -log(p);
		score(bump(ii,jj):30,:) = smooth_2D(score(bump(ii,jj):30,:),[15, 1, 2], [2, 0.2, 2]);

		%shink
		score = score(1:dis_cut,:); %30x16 -> 24x16
		score(score==0) = 0; %fix -0
		score(score>range) = range;
		score(score<-range) = -range;

		%debug
		if debug
			printf("score(%d,%d): max=%f, min=%f\n", ii, jj, max(max(score)), min(min(score)));
			contourf(xx,yy,score);
			fname = [ "score_", int2str(ii), "_", int2str(jj), ".png" ];
			print(fname, "-dpng")
		else
			fprintf(fp_ang, "ANGLE1: %s %s -0.25 12.25 26 -1.0625 1.0625 18\n", list(ii,:), list(jj,:));
			for yy=1:18
				fprintf(fp_ang, "%.6f ", 0);
			end
			fprintf(fp_ang, "\n");
			for xx=1:dis_cut
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
	end
end

if debug
	printf("Done!\n");
	set(f, "visible", "on")
else
	fclose(fp_ang);
end

