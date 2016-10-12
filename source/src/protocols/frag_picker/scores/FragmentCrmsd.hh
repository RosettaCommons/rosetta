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
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_FragmentCrmsd_hh
#define INCLUDED_protocols_frag_picker_scores_FragmentCrmsd_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// mini

#include <core/types.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

//Auto Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <iostream>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<core::Real> > Matrix;

/// @brief  scores a fragment by its crmsd to the given reference structure
class FragmentCrmsd: public CachingScoringMethod {
public:

	/// @brief  creates a crmsd-based scoring function.
	/// @details fragments will be compared to a given pose, which should have the same number of residues a the query sequence
	FragmentCrmsd(core::Size, core::Real, bool, core::pose::PoseOP);

	/// @brief  creates a crmsd-based scoring function.
	/// @details fragments will be compared to given coordinates, which should have the same number of residues a the query sequence
	FragmentCrmsd(core::Size, core::Real, bool, utility::vector1< utility::vector1<core::Real> >);

	~FragmentCrmsd();

	void do_caching(VallChunkOP);
	void clean_up();
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	core::Size n_atoms_;
	core::pose::PoseOP reference_pose_;
	ObjexxFCL::FArray2D_double reference_coordinates_;
	ObjexxFCL::FArray2D_double chunk_coordinates_;
	ObjexxFCL::FArray2D_double fragment_coordinates_;
	ObjexxFCL::FArray2D_double frag_pos_coordinates_;
	ObjexxFCL::FArray1D_double weights_;

	void fill_CA_coords(core::pose::Pose const &, ObjexxFCL::FArray2_double &, core::Size);
};

/// @brief  Maker class that produces a new FragmentCrmsd object
class MakeFragmentCrmsd: public MakeFragmentScoringMethod {
public:

	MakeFragmentCrmsd() :
		MakeFragmentScoringMethod("FragmentCrmsd") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_FragmentCrmsd_HH
