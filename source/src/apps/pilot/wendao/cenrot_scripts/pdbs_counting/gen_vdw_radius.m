[nonidcen_counts] = load_nonidcen_table();


%cutoff = 0.003;
%pcut = 0.1

aaa = 16;
cutoff = 0.01;
pcut = 1;

%clean the clash
nonidcen_counts(:,:,1:6) = 0;

dic = [ "Nbb";
		"CAbb";
		"CB";
		"CObb";
		"OCbb";
		"ALA";
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


rvdw1 = 0;
rvdw2 = 0;
for aandx = aaa
printf("%s ",dic(aandx,:));
for ii=1:25
	subplot(5,5,ii);
	%r = 0.25:0.5:15;
	r = 0.1:0.2:10;
	shell = (4.0*pi*r.*r)';
	Nr = squeeze(nonidcen_counts(ii, aandx, :));
	Nint = Nr;
	Pr = Nr ./ shell;
	Nall = 0;
	for ndx = 1:50
		Nall = Nall+Nr(ndx);
		Nint(ndx) = Nall;
	end
	Nc = Nall*cutoff;

	rm = 0:0.01:10;
	Nrm = interp1( r, Nint, rm, 'linear');
	Prm = interp1( r, Pr, rm, 'linear');

	oldp = 0;
	P1=-1;
	P2=-1;
	for ndx=2:10000
		if Nrm(ndx) > Nc
			rvdw1 = rm(ndx);
			break;
		end
		if Prm(ndx) > pcut
			rvdw2 = rm(ndx);
			break;
		end
	end
	rvdw = rvdw1;
	if rvdw<rvdw2
		rvdw = rvdw2;
	end

	printf("%.3f ",rvdw*rvdw);
	rx = 0:0.1:10;
	Px = interp1( r, Pr, rx, 'linear');

	Pvdw = exp(-(rvdw**2-rx.*rx).**2/rvdw**2);
	Pvdw(rx>rvdw) = 1;

	%%draw the g(r) picture
	h=plot( rm, Prm, 'r', rx, Pvdw, 'b');
	axis([0 10 0 1.5]);
	t = dic(ii,:);
	text(0, Nall/100000, t);
	set (h(1), "linewidth", 2)
	set (h(2), "linewidth", 2)
end
printf("\n");
end

