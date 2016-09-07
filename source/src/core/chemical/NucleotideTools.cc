// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/NucleotideTools.cc
/// @author Christoffer Norn (ch.norn@gmail.com)

#include <core/chemical/NucleotideTools.hh>
#include <map>
#include <core/types.hh>
#include <string>
#include <utility/vector0.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>


namespace core {
namespace chemical {
namespace NucleotideTools {

using namespace std;

string
codon2aa( string const & codon ) {
	std::string codons[64] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG", "CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG",
		"GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"};
	string aminoAcids[64]={"F","F","L","L", "S","S","S","S","Y","Y", "Stop","Stop", "C","C","Stop", "W",
		"L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R",
		"I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R",
		"V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"};
	map< string, string> nt_2_aa;
	for ( int i=0; i<64; i++ ) {
		nt_2_aa[codons[i]] = aminoAcids[i];
	}
	return nt_2_aa[codon];
}

string
aa2randomCodon( char const & aa ) {
	string codons[61] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG",
		"CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT",
		"CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC",
		"AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA",
		"GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"};
	char aminoAcids[61] = {'F','F','L','L','S','S','S','S','Y','Y','C','C','W','L','L','L','L','P','P','P',
		'P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K',
		'K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G',
		'G'};

	// Setup aa to codon table map
	map< char, vector<string> > aa_2_nt;
	for ( int i=0; i<61; i++ ) {
		aa_2_nt[ aminoAcids[i] ].push_back( codons[i] );
	}

	vector<string> nts = aa_2_nt[ aa ];

	numeric::random::random_permutation(nts, numeric::random::rg());
	return nts[0];
}


} // NucleotideTools
} // chemical
} // core
