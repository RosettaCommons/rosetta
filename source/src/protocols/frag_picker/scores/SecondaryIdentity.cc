// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/SecondaryIdentity.cc
/// @brief  scores a fragment by an amino acid sequence identity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/SecondaryIdentity.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace ObjexxFCL::format;

bool SecondaryIdentity::score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {

	core::Real totalScore = 0;
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		VallResidueOP r = f->get_residue(i);
		if ( r->ss() == query_[i + f->get_first_index_in_query() - 2] ) {
			totalScore += REWARD;
		} else {
			totalScore += PENALTY;
		}
	}
	totalScore /= (core::Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

bool SecondaryIdentity::describe_score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map, std::ostream& out) {

	core::Real totalScore = 0;

	out << f->get_chunk()->get_pdb_id() << " " << I(5,
		f->get_first_index_in_vall()) << " ";
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		out << f->get_residue(i)->ss();
	}
	out << std::endl << "            ";
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		VallResidueOP r = f->get_residue(i);
		if ( r->ss() == query_[i + f->get_first_index_in_query() - 2] ) {
			totalScore += REWARD;
			out << "+";
		} else {
			totalScore += PENALTY;
			out << "-";
		}
	}
	out << "\nquery " << I(5, f->get_first_index_in_query()) << " ";
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		out << query_[i + f->get_first_index_in_query() - 2];
	}
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		out << "\nTotal score " << F(5, 3, totalScore) << " ACCEPTED"
			<< std::endl;
	} else {
		out << "\nTotal score " << F(5, 3, totalScore) << " REJECTED"
			<< std::endl;
	}
	totalScore /= (core::Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

} // scores
} // frag_picker
} // protocols
