cutoff = 12;

%load all table data
[dic, pair_table, env_table] = load_pair_env_table();
[na1, na2, nb] = size(pair_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gen pair score file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_pair = fopen("cen_rot_pair.log", "w");
for ii=1:na1
	for jj=ii:na2
		fprintf(fp_pair, "PAIR: %s %s ", dic(ii,:), dic(jj,:));
		fprintf(fp_pair, "%.2f %.2f ", 0.25, 0.25+(cutoff*2-1)*0.5);
		[x,y1,y2]=calc_cenpair_scores(pair_table, ii, jj, cutoff);

		fprintf(fp_pair, "%d ", length(y2));
		for k=1:length(y2)
			fprintf(fp_pair, "%f ", y2(k));
		end
		fprintf(fp_pair, "\n");
	end
end
fprintf(fp_pair, "# END\n");
fclose(fp_pair);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gen env score file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%smooth data
env_counts_smoothed=zeros(20,40);
L=4; %half window
CONV1 = gausswin(2*L+1, 2);
CONV1 = CONV1/sum(CONV1);
for i = 1:20
    smoothed = conv(squeeze(env_table(i,:)), CONV1);
    env_counts_smoothed(i,:) = smoothed(L+1:end-L);
end

env_sums_aa = sum(env_counts_smoothed, 2);
%P_env=zeros(20,40);
%for i = 1:40
%    P_env(:,i) = env_counts_smoothed(:,i) ./ env_sums_aa;
%end
%env_sums_tot = squeeze(sum(P_env));
%for i = 1:20
%    P_env(i,:) = P_env(i,:) ./ env_sums_tot * 20;
%end
%score = -log(P_env);

sums_env = sum(env_counts_smoothed); %1x40
P_env = zeros(20, 40);
for i = 1:20
    P_env(i,:) = env_counts_smoothed(i,:) ./ sums_env;
end
totalaa = sum(sums_env);
P_aa = env_sums_aa / totalaa;
for i = 1:20
    score(i,:) = -log(P_env(i,:))+log(P_aa(i));
end

fp_env = fopen("cen_rot_env.log", "w");
for ii=1:20
    x = [1:40];
    y1 = score(ii,x);
    y1(1:5) = y1(6);
	if (dic(ii,:)=="ARG" || dic(ii,:)=="LYS" || dic(ii,:)=="GLU")
		y1(1:8) = y1(9);
	endif
    y1(35:40) = y1(34);

    %output
    fprintf(fp_env, "ENV: %s ", dic(ii, :));
    for jj=1:40
        fprintf(fp_env, "%.6f ", y1(jj));
    end
    fprintf(fp_env, "\n");
end
fprintf(fp_env, "# END");
fclose(fp_env);

