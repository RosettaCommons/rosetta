env_counts = load_env_table();
env_ref = load_env_table_ref();

%smooth data
env_counts_smoothed=zeros(20,40);
env_counts_smoothed_ref=zeros(20,40);
L=6; %half window
CONV1 = gausswin(2*L+1, 3);
CONV1 = CONV1/sum(CONV1);
for i = 1:20
    smoothed = conv(squeeze(env_counts(i,:)), CONV1);
    smoothed_ref = conv(squeeze(env_ref(i,:)), CONV1);
    env_counts_smoothed(i,:) = smoothed(L+1:end-L);
    env_counts_smoothed_ref(i,:) = smoothed_ref(L+1:end-L);
end

env_sums_aa = sum(env_counts_smoothed, 2);
env_sums_aa_ref = sum(env_counts_smoothed_ref, 2);

P_env=zeros(20,40);
P_env_ref=zeros(20,40);
for i = 1:40
	P_env(:,i) = env_counts_smoothed(:,i) ./ env_sums_aa;
	P_env_ref(:,i) = env_counts_smoothed_ref(:,i) ./ env_sums_aa_ref;
end
env_sums_tot = squeeze(sum(P_env));
env_sums_tot_ref = squeeze(sum(P_env_ref));
for i = 1:20
	P_env(i,:) = P_env(i,:) ./ env_sums_tot * 20;
	P_env_ref(i,:) = P_env_ref(i,:) ./ env_sums_tot_ref * 20;
end
score = -log(P_env);
score_ref = -log(P_env_ref);

dic = [ "GLU";
	"LYS";
	"SER";
	"ILE";
	"ASN";
	"CYS";
	"ASP";
	"LEU";
	"GLY";
	"VAL";
	"ARG";
	"PHE";
	"PRO";
	"THR";
	"TYR";
	"ALA";
	"MET";
	"HIS";
	"GLN";
	"TRP" ];

dic2 = [ "ALA";
        "CYS";
        "ASP";
        "GLU";
        "PHE";
        "GLY";
        "HIS";
        "ILE";
        "LYS";
        "LEU";
        "MET";
        "ASN";
        "PRO";
        "GLN";
        "ARG";
        "SER";
        "THR";
        "VAL";
        "TRP";
        "TYR" ];

[cbeta6,cbeta12,cenpack,env,pair_counts,pair_old] = load_tables();

for ii=1:20
    subplot(4,5,ii);
	x = [1:40];
	y1 = score(ii,x);
	y1_ref = score_ref(ii,x);
	y1(1:3) = y1(4);
	y1_ref(1:3) = y1_ref(4);
	y1(36:40) = y1(36);
	y1_ref(36:40) = y1_ref(36);
	for jj=1:20
		if (dic2(jj,:)==dic(ii,:))
			break;
		end
	end
    h=plot(x, env(jj,x), 'b', x, y1_ref, 'g', x, y1, 'r');
    t = dic(ii,:);
    text(0, y1(1), t);
    set (h(1), "linewidth", 2)
    set (h(2), "linewidth", 2)
    set (h(3), "linewidth", 2)
	set (gca, 'xtick', [0:10:30])
	if (ii>15)
		set (gca, 'xticklabel', {'0','10','20','30'}) 
	else
		set (gca, 'xticklabel', {}) 
	end
	set (gca, 'fontsize', 4) 

	%output
	printf("ENV_LOG: %s ", dic(ii, :));
	for jj=1:40
		printf("%.6f ", y1(jj));
	end
	printf("\n");
end
