// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CustomScoreSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/CustomScoreSelector.hh>

#include <protocols/frag_picker/FragmentSelectingRule.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <utility/vector1.hh>
// utility headers
#include <core/types.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace frag_picker {

static THREAD_LOCAL basic::Tracer trCustomScoreSelector( "protocols.frag_picker.CustomScoreSelector" );

typedef utility::vector1<std::pair< ScoredCandidate, core::Real> > ExtraScoreVector;

bool sort_function( ExtraScoreVector::value_type p1, ExtraScoreVector::value_type p2 ) {
	return ( p1.second < p2.second );
}

void CustomScoreSelector::select_fragments(
	ScoredCandidatesVector1 const& input_canditates,
	ScoredCandidatesVector1& output_selection )
{
	Size n = frags_per_pos();
	if ( n > input_canditates.size() ) {
		n = input_canditates.size();
	}
	trCustomScoreSelector.Debug << "Selecting " << n << "fragments from "
		<< input_canditates.size() << " candidates" << std::endl;
	ExtraScoreVector tmp;
	scores::FragmentScoreMapOP map = scoring_scheme_->create_empty_map();
	for ( Size i=1; i<=input_canditates.size(); i++ ) {
		scoring_scheme_->do_caching(input_canditates[i].first->get_chunk());
		scoring_scheme_->score_fragment( input_canditates[i].first, map );
		core::Real v = scoring_scheme_->total_score( map );
		scoring_scheme_->clean_up();
		ExtraScoreVector::value_type p(input_canditates[i], v);
		tmp.push_back( p );
	}

	std::sort(tmp.begin(), tmp.end(), sort_function);
	for ( Size i = 1; i <= n; i++ ) {
		output_selection.push_back( tmp[i].first );
	}
}


} // frag_picker
} // protocols

