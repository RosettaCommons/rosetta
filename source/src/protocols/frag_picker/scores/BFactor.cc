// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/BFactor.cc
/// @brief  a base class for fragment scoring
/// @author rvernon@u.washington.edu

#include <protocols/frag_picker/scores/BFactor.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

bool BFactor::score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {

	// describe_score(f, empty_map, std::cerr);

	Real totalScore = 0;
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		VallResidueOP r = f->get_residue(i);
		totalScore += r->bF();
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

bool BFactor::describe_score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map, std::ostream& out) {

	Real totalScore = 0;

	out << f->get_chunk()->get_pdb_id() << "  " << I(5,
		f->get_first_index_in_vall()) << " ";
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		out << f->get_residue(i)->bF();
	}
	out << std::endl;

	//THIS FUNCTION IS POINTLESS!!!

	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

} // scores
} // frag_picker
} // protocols


