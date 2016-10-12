// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.hh
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Nikolas Sgourakis sgourn@u.w.edu

#ifndef INCLUDED_protocols_frag_picker_scores_JCoupling_hh
#define INCLUDED_protocols_frag_picker_scores_JCoupling_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// mini
#include <core/types.hh>


#include <protocols/frag_picker/JCouplingIO.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {


/// @brief  scores a fragment by the JCouplings
class JCoupling: public CachingScoringMethod {
public:

	/// @brief  creates a
	/// @details
	JCoupling(core::Size, core::Real, bool, JCouplingIO&);

	void do_caching(VallChunkOP);
	void clean_up();
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);
	bool score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	//std::string cached_scores_id_; // cache is not yet built
	JCouplingIO data_;
	core::Size len_;
	core::Real A_, B_, C_, THETA_;
};

/// @brief  Matker class that produces a new JCoupling object
class MakeJCoupling: public MakeFragmentScoringMethod {
public:

	MakeJCoupling() :
		MakeFragmentScoringMethod("JCoupling") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_JCoupling_HH
