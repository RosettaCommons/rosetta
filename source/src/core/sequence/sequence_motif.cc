// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sequence/util.cc
/// @brief Utilities for working with sequence motifs.
/// @author Jared Adolf-Bryfogle

// Unit header
#include <core/sequence/sequence_motif.hh>

// C/C++ headers
#include <algorithm>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <utility/string_util.hh>
#include <string>

namespace core {
namespace sequence {

using namespace utility;

static basic::Tracer tr( "core.sequence.sequence_motif" );

utility::vector1<std::string>
split_sequence_motif( std::string const & motif ){

	//This is because NV means something different than [NV].  This enables us to detect that and account for it.
	// SO - any string that has + in it does not come from [NV].

	std::string s = replace_in( motif, "[", "+[" );
	s = replace_in( s, "]", "]+");

	s = strip( s, ']');
	s = strip( s, '[');

	s = replace_in( s, "][", ",");
	s = replace_in( s, "]", ",");
	s = replace_in( s, "[", ",");

	utility::vector1< std::string > fsplit = string_split( s, ',');
	utility::vector1< std::string > ssplit;
	for ( std::string str : fsplit ) {
		if ( ( str != "++") && str != ("+") ) {
			ssplit.push_back(str);
		}
	}

	//If nothing needs splitting, it is per-position.
	if ( ssplit.size() == 1 ) {
		ssplit[1] = "+"+ssplit[1]+"+";
	}

	vector1< std::string > final_split_motif;
	for ( std::string str : ssplit ) {
		if ( contains(str, "+") ) {
			std::string rem_str = remove_from_string(str, "+");
			for ( char c : rem_str ) {
				final_split_motif.push_back( to_string(c) );
			}
		} else {
			final_split_motif.push_back(str);
		}
	}
	return final_split_motif;
}

core::Size
get_motif_length( std::string const & motif ){
	vector1< std::string > final_split_motif = split_sequence_motif(motif);
	return final_split_motif.size();
}

std::string
get_design_sequence_motif_syntax(){

	std::string motif_str =
		"  This is slightly similar to a regex, but not quite. We are not matching a sequence,\n"
		"   we are designing in a motif regardless of the current sequence, anywhere in a protein.\n"
		"\n"
		"   - Each letter corresponds to a position. Using [ ] indicates a more complicated expression for that position.\n"
		"   - An X indicates it can be anything, and that we are designing here.\n"
		"   - An AA Letter, like V, indicates that that position will be designed to a V.\n"
		"   - A - charactor indicates that that position stays with whatever it is currently.  We essentially skip this position.\n"
		"   - An expression like: [^PAV] indicates that we will design anything except Proline, Alanine, and Valine \n"
		"   - An expression like: [NTS] indicates that that position can be Asparigine, Threonine, or Serine and \n"
		"      only of these will be enabled during the design.\n"
		"   - RESFILE commands are accepted as well. These require a % charactor in from of the whole expression.\n"
		"     For example [%POLAR] would set that position to only polar design.\n"
		"     This is exactly the same as a resfile line, so you can even do NC like so: \n"
		"      [%EMPTY NC R2 NC T6 NC OP5]\n"
		"\n"
		" EXAMPLE:\n"
		"  Glycosylation N-Linked motif design: N[^P][ST]\n"
		"\n";

	return motif_str;
}



} // sequence
} // core
