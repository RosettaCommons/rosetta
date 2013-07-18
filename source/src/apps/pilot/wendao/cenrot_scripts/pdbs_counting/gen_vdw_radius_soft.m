[nonidcen_counts] = load_nonidcen_table();

%clean the clash
%nonidcen_counts(:,:,1:6) = 0;

dic = [ "Nbb";
		"CAbb";
		"CB";
		"CObb";
		"OCbb";
		"CEN_ALA";
        "CEN_CYS";
        "CEN_ASP";
        "CEN_GLU";
        "CEN_PHE";
        "CEN_GLY";
        "CEN_HIS";
        "CEN_ILE";
        "CEN_LYS";
        "CEN_LEU";
        "CEN_MET";
        "CEN_ASN";
        "CEN_PRO";
        "CEN_GLN";
        "CEN_ARG";
        "CEN_SER";
        "CEN_THR";
        "CEN_VAL";
        "CEN_TRP";
        "CEN_TYR" ];

cutoff = 0.001;

for aandx = aaa
printf("%s ",dic(aandx,:));
for ii=1:25
	subplot(5,5,ii);
	%r = 0.25:0.5:15;
	r = 0.1:0.2:10;
	shell = (4.0*pi*r.*r)';
	Nr = squeeze(nonidcen_counts(ii, aandx, :));
	Nint = Nr;
	Pr = Nr ./ shell / 0.2;
	Nall = 0;
	for ndx = 1:50
		Nall = Nall+Nr(ndx);
		Nint(ndx) = Nall;
	end
	Nc = Nall*cutoff;

	rm = 0:0.001:10;
	sm = 4.0*pi*rm.*rm;
	Nrm = interp1( r, Nint, rm, 'linear');
	Prm = interp1( r, Pr, rm, 'linear');

	p_end = Pr(50)/3;
	for ndx=2:10000
		if Prm(ndx)>p_end
			break
		end
	end
	rvdw = rm(ndx);

	printf("%.3f ",rvdw*rvdw);
	rx = 0:0.1:10;
	Px = interp1( r, Pr, rx, 'linear');

	Pvdw = exp(-(rvdw**2-rx.*rx).**2/rvdw);
	Pvdw(rx>rvdw) = 1;

	%%draw the g(r) picture
	h=plot( rx, Px, 'r', rx, Pvdw*p_end, 'b');
	t = dic(ii,:);
	text(0, Nall/100000, t);
	set (h(1), "linewidth", 2)
	set (h(2), "linewidth", 2)
end
printf("\n");
end

