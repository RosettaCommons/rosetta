// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/CachingScoringMethod.hh
/// @brief  adds cashing functionality to a FragmentScoringMethod object
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_CachingScoringMethod_hh
#define INCLUDED_protocols_frag_picker_scores_CachingScoringMethod_hh

// type headers
#include <core/types.hh>
// package headers
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <protocols/frag_picker/VallChunk.fwd.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <utility/exit.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// APL: Always always declare your OP typedefs if you're declaring a polymorphic class.
/// I shouldn't have to do this for you.
class CachingScoringMethod;
typedef utility::pointer::shared_ptr< CachingScoringMethod > CachingScoringMethodOP;
typedef utility::pointer::shared_ptr< CachingScoringMethod const > CachingScoringMethodCOP;


/// @brief
class CachingScoringMethod: public FragmentScoringMethod {
public:

	CachingScoringMethod(Size priority, Real lowest_acceptable_value,
			bool use_lowest, std::string name) :
		FragmentScoringMethod(priority, lowest_acceptable_value, use_lowest, name) {
	}

	virtual void do_caching(VallChunkOP) = 0;
	virtual void clean_up() = 0;
	virtual bool cached_score(FragmentCandidateOP, FragmentScoreMapOP) {
	  utility_exit_with_message( "ERROR: unimplemented cached_score() method. Your score function could not use cache\n" );
	  return true;
	}
	virtual bool score(FragmentCandidateOP fragment,
	                FragmentScoreMapOP scores) {
			return score(fragment,scores); };
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_CachingScoringMethod_HH */
