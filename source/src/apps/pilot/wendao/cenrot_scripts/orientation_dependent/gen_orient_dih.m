[list, dat_ang1, dat_ang2, dat_dih] = load_angle_r_table();

debug = 0;
cutoff = 10; %bin tol
dis_cut = 24; %pair term cutoff 12A
range = 0.6;
ntype = size(list)(1);

dat_all = dat_dih;
dat_a_a_r_dih = zeros(ntype, ntype, 30, 18);
bump = zeros(ntype, ntype);
for aa = 1:ntype
    for bb = 1:ntype
        dat_all(aa,bb,:,:) += dat_dih(bb,aa,:,:);
        tmpdata = squeeze(dat_all(aa,bb,:,:));
        tmpsumr = sum(tmpdata, 2);
        for cc = 1:30
            if tmpsumr(cc)>cutoff*18
                bump(aa,bb) = cc;
                break;
            end
        end
        dat_a_a_r_dih(aa,bb,cc:30,:) = dat_all(aa,bb,cc:30,:);
    end
end

sum_a_a_r = sum(dat_a_a_r_dih,4); %sum over ang
dat_a_r_ang1 = squeeze(sum(dat_a_a_r_dih, 2)); %sum over the second aa
sum_a1_r = sum(dat_a_r_ang1, 3); %sum over ang
dat_a_r_ang2 = squeeze(sum(dat_a_a_r_dih)); %sum over the first aa
sum_a2_r = sum(dat_a_r_ang2, 3); %sum over ang

dat_r_ang = squeeze(sum(dat_a_r_ang1)); %sum over all aa
sum_r = sum(dat_r_ang,2); %sum over ang

p_a_a_r_ang = zeros(ntype, ntype, 30, 18);
p_a_r_ang1 = zeros(ntype, 30, 18);
p_a_r_ang2 = zeros(ntype, 30, 18);
p_r_ang = zeros(30, 18);

for kk = 1:30
    for ii = 1:ntype
        for jj = 1:ntype
            if sum_a_a_r(ii,jj,kk)>cutoff*18
                p_a_a_r_ang(ii,jj,kk,:) = dat_a_a_r_dih(ii,jj,kk,:) / sum_a_a_r(ii,jj,kk);
            else
                p_a_a_r_ang(ii,jj,kk,:) = 1.0/18;
            end
        end
        if sum_a1_r(ii,kk)>cutoff*18*ntype
            p_a_r_ang1(ii,kk,:) = dat_a_r_ang1(ii,kk,:) / sum_a1_r(ii,kk);
        else
            p_a_r_ang1(ii,kk,:) = 1.0/18;
        end
        if sum_a2_r(ii,kk)>cutoff*18*ntype
            p_a_r_ang2(ii,kk,:) = dat_a_r_ang2(ii,kk,:) / sum_a2_r(ii,kk);
        else
            p_a_r_ang2(ii,kk,:) = 1.0/18;
        end
    end
    if sum_r(kk)>cutoff*18*ntype*ntype
        p_r_ang(kk,:) = dat_r_ang(kk,:) / sum_r(kk);
    else
        p_r_ang(kk,:) = 1.0/18;
    end
end

p_r_ang = smooth_2D(p_r_ang, [15, 1, 2], [360, 20, 2]);

%make up a prior count to force the score down to 0
prior = zeros(30,18);
for i = 15:30
    prior(i,:) = (i-15)*0.02/15;
end

%debug
if debug
    [xx, yy] = meshgrid(1:18, 1:24);
    f = figure();
    set(f, "visible", "off")
else
    fp_ang = fopen("cen_rot_pair_dih.log", "w");
end

for ii = 1:ntype
    for jj = ii:ntype
        %skip GLY/ALA
        if bump(ii,jj)==0
            continue;
        end

        mask0 = squeeze(dat_all(ii,jj,:,:))<cutoff; %no enough data
        p_ii_jj = ones(30, 18)/18;
        p_ii_jj(bump(ii,jj):30,:) = squeeze(p_a_a_r_ang(ii, jj, bump(ii,jj):30, :));
        mask1 = p_ii_jj==0;
        p_ii_jj(bump(ii,jj):30,:) = smooth_2D(p_ii_jj(bump(ii,jj):30,:),[15, 1, 2], [360, 20, 1]);
        p_ii = ones(30, 18)/18;
        p_jj = ones(30, 18)/18;
        p_ii(bump(ii,jj):30,:) = p_a_r_ang1(ii,bump(ii,jj):30,:);
        p_jj(bump(ii,jj):30,:) = p_a_r_ang2(jj,bump(ii,jj):30,:);
        p_ii_p_jj = p_ii .* p_jj;
        mask2 = p_ii_p_jj==0;
        p_ii_p_jj(bump(ii,jj):30,:) = smooth_2D(p_ii_p_jj(bump(ii,jj):30,:),[15, 1, 2], [360, 20, 1]);

        mask = (mask0+mask1+mask2)>0;
        p_ii_jj(mask)=1.0/18;
        p_ii_p_jj(mask)=1.0/18/18;

        p = (p_ii_jj .* p_r_ang + prior) ./ (p_ii_p_jj + prior);
        p(mask) = 1;
        p(1:bump(ii,jj)-1,:) = 1;
        score = -log(p);
        score(bump(ii,jj):30,:) = smooth_2D(score(bump(ii,jj):30,:),[15, 1, 2], [360, 20, 1]);

        %shink
        score = score(1:dis_cut,:); %30x18 -> 24x18
        score(score==0) = 0; %fix -0
        score(score>range) = range;
        score(score<-range) = -range;

        %debug
        if debug
            printf("score(%d,%d): max=%f, min=%f\n", ii, jj, max(max(score)), min(min(score)));
            contourf(xx,yy,score);
            fname = [ "dih_score_", int2str(ii), "_", int2str(jj), ".png" ];
            print(fname, "-dpng")
        else
            fprintf(fp_ang, "DIHEDRAL: %s %s 0.25 11.75 24 -170 170 18\n", list(ii,:), list(jj,:));
            for xx=1:dis_cut
                for yy=1:18
                    fprintf(fp_ang, "%.6f ", score(xx,yy));
                end
                fprintf(fp_ang, "\n");
            end
        end
    end
end

if debug
    printf("Done!\n");
    set(f, "visible", "on")
else
    fclose(fp_ang);
end

