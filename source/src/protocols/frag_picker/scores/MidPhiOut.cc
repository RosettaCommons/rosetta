// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.cc
/// @brief  Object that scores a fragment by root mean square deviation of Phi && Psi dihedrals
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/MidPhiOut.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// option key includes
#include <basic/options/keys/OptionKeys.hh>


// utils
#include <ObjexxFCL/FArray1D.hh>

// C++

// Boost

#include <utility/vector1.hh>

#ifdef WIN32
#include <protocols/frag_picker/FragmentPicker.hh>
#endif


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

MidPhiOut::MidPhiOut(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "MidPhiOut") {
}

void MidPhiOut::do_caching(VallChunkOP current_chunk) {

	chunk_phi_.redimension(current_chunk->size());
	for ( core::Size i = 1; i <= current_chunk->size(); ++i ) {
		VallResidueOP r = current_chunk->at(i);
		chunk_phi_(i) = r->phi();
	}
}

bool MidPhiOut::score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool MidPhiOut::cached_score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {

	std::string tmp = fragment->get_chunk()->chunk_key();

	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(fragment->get_chunk());
		cached_scores_id_ = tmp;
	}

	core::Size offset_v = fragment->get_first_index_in_vall() - 1;

	auto r_phi = static_cast< core::Real > ( fragment->get_length() );
	r_phi = (r_phi/2)+0.5+offset_v;
	auto i_phi = static_cast< core::Size > ( r_phi );

	core::Real phi = chunk_phi_(i_phi);

	scores->set_score_component( phi, id_);

	return true;
}

void MidPhiOut::clean_up() {
}

/// @brief Creates a MidPhiOut scoring method
/// @param priority - priority of the scoring method. The higher value the earlier the score
///  will be evaluated
/// @param lowest_acceptable_value - if a calculated score is higher than this value,
///  fragment will be neglected
/// @param FragmentPickerOP object - not used
/// @param line - the relevant line extracted from the scoring configuration file that defines this scoring method
///   It could look like: "MidPhiOut                140     -5.0     100.0 additional_string"
///  where 140, -5.0 && 100.0 are priority, weight && treshold, respectively.
///  The additional string may be:
///  - empty: then the maker tries to create a scoring object from a TALOS file
///   trying in::file::talos_phi_psi flag. If fails, will try to use a pose from in::file::s
///  - a pdb file, pdb extension is necessary. This will create a pose && steal Phi && Psi
///  - a TALOS file with Phi/Psi prediction (tab extension is necessary)
FragmentScoringMethodOP MakeMidPhiOut::make(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
	, std::string ) {

	return (FragmentScoringMethodOP) FragmentScoringMethodOP( new MidPhiOut(priority,
		lowest_acceptable_value, use_lowest) );
}

} // scores
} // frag_picker
} // protocols
