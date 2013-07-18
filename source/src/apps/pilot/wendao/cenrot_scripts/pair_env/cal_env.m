[dic, pairtab, env_counts] = load_pair_env_table();

%smooth data
env_counts_smoothed=zeros(20,40);
L=4; %half window
CONV1 = gausswin(2*L+1, 2);
CONV1 = CONV1/sum(CONV1);
for i = 1:20
    smoothed = conv(squeeze(env_counts(i,:)), CONV1);
    env_counts_smoothed(i,:) = smoothed(L+1:end-L);
end

%P(env|aa)/P(env)
env_sums_aa = sum(env_counts_smoothed, 2);
P_env=zeros(20,40);
for i = 1:40
	P_env(:,i) = env_counts_smoothed(:,i) ./ env_sums_aa;
end
env_sums_tot = squeeze(sum(P_env));
for i = 1:20
	P_env(i,:) = P_env(i,:) ./ env_sums_tot * 20;
end
score1 = -log(P_env);

%P(aa|env)/P(aa)
sums_env = sum(env_counts_smoothed); %1x40
P_env = zeros(20, 40);
for i = 1:20
	P_env(i,:) = env_counts_smoothed(i,:) ./ sums_env;
end
totalaa = sum(sums_env);
P_aa = env_sums_aa / totalaa;
for i = 1:20
	score2(i,:) = -log(P_env(i,:))+log(P_aa(i));
end

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
	y1 = score1(ii,x);
	y1(1:5) = y1(6);
	y1(35:40) = y1(34);
	y2 = score2(ii,x);
	y2(1:5) = y2(6);
	y2(35:40) = y2(34);
	for jj=1:20
		if (dic2(jj,:)==dic(ii,:))
			break;
		end
	end
    h=plot(x, env(jj,x), 'g', x, y1, 'r', x, y2, 'b');
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
		printf("%.6f ", y2(jj));
	end
	printf("\n");
end
