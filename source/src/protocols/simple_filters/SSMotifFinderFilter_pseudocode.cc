Pose template_pose
Size template_stem1,2
Size from_res, to_res
string motif
Real rmsd
string filename
jump

the algorithm:

{  
std::string const SS_sequence( dssp( pose ) ); //this is how I use dssp? pose is the structure loaded from the disc?
[vector of strings? what's the output of regex_search?] motif_matches = regex_search (SS_sequence.begin(), SS_sequence.end(), match_results, motif, flag); //[apply 'regex_search' on 'SS_sequence' to search for 'motif' regex]. what is match_results, flag?
for [loop over each string 'match' in 'motif_matches']{
	for (Size i=[match_begining]+from_res ; i<= [match_begining]+to_res; i++){  
		if (to_res > Size [match_end])
        		to_res = [match_end];
		int counter=0;
		counter++;
		[create array 'matches_splice_sites' of Size** in the size of [2,counter*(amount of matches (char*?) in motif_matches)];
		save [match_begining,i] in the array;  //make sure to save the chain data
	};
};
[delete repeating pairs of [match_begining,i] in 'matches_splice_sites' so get each point once. maybe through the 'sort' option];
for [go through each pair in the array 'match_splice_sites']{
	Residue s2**=[copy the second residue in the pair];  	//how is the Residue variable called?
	apply 'jump' from [residue in place (1st residue in the pair) to residue s2**];
	float rmsd_apparent = [calculate rmsd between residues S2** and [2nd residue in the pair]];
	if ('rmsd_apparent' <= 'rmsd')
		save [ [PDB name], [1st residue in the pair], [2nd residue in the pair], 'rmsd_apparent' \n] in 'filename';
};
};

