[cbeta6,cbeta12,cenpack,env,pair_counts,pair_old] = load_tables();
[dic2, pairtab, envtab] = load_pair_env_table();
%newcen_counts = squeeze(pairtab(1,:,:,:));
newcen_counts = pairtab;
cutoff=12;
dic = [ "ALA";
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

aaii = 0;
for jj=1:20
	if (dic2(jj,:)==dic(aandx,:)) 
		aaii=jj;
		break;
	endif
end

for ii=1:20
	subplot(4,5,ii);

	for jj=1:20
		if (dic2(jj,:)==dic(ii,:)) 
			break;
		endif
	end

	%[x,y1,y2]=calc_cenpair_scores(pair_counts, aandx, ii, cutoff);
	[x,y3,y4]=calc_cenpair_scores(newcen_counts, aaii, jj, cutoff);
	yold = calc_ros_pair(pair_old, 15, ii);
	%h=plot( x, y1, 'r', x, y3, 'g', 3:.1:12, yold, 'b');
	h=plot( x, y3, 'g', 3:.1:12, yold, 'b', x, y4, 'r');

	axis([4 12 0 1], "autoy"); 
	t = dic(ii,:);
	text(4, y3(1), t);
	set (h(1), "linewidth", 2)
	set (h(2), "linewidth", 2)
	set (h(3), "linewidth", 2)
end

