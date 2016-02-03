// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/FragmentCrmsdResDepth.hh
/// @brief  Object that scores a fragment by its crmsd and residue depth to the native
/// @author David E Kim

#ifndef INCLUDED_protocols_frag_picker_scores_FragmentCrmsdResDepth_hh
#define INCLUDED_protocols_frag_picker_scores_FragmentCrmsdResDepth_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

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

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  scores a fragment by its crmsd to the given reference structure
class FragmentCrmsdResDepth: public CachingScoringMethod {
public:

	/// @brief  creates a crmsd-based scoring function.
	/// @details fragments will be compared to a given pose, which should have the same number of residues a the query sequence
	FragmentCrmsdResDepth(Size, Real, bool, core::pose::PoseOP, utility::vector1<core::Real> & query_residue_depth);

	/// @brief  creates a crmsd-based scoring function.
	/// @details fragments will be compared to given coordinates, which should have the same number of residues a the query sequence
	FragmentCrmsdResDepth(Size, Real, bool, utility::vector1< utility::vector1<Real> >,  utility::vector1<core::Real> & query_residue_depth);

	~FragmentCrmsdResDepth();

	void do_caching(VallChunkOP);
	void clean_up();
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	Size n_atoms_;
	pose::PoseOP reference_pose_;
	ObjexxFCL::FArray2D_double reference_coordinates_;
	ObjexxFCL::FArray2D_double chunk_coordinates_;
	ObjexxFCL::FArray2D_double fragment_coordinates_;
	ObjexxFCL::FArray2D_double frag_pos_coordinates_;
	ObjexxFCL::FArray1D_double weights_;
	utility::vector1<core::Real> query_residue_depth_;

	void fill_CA_coords(core::pose::Pose const &, ObjexxFCL::FArray2_double &, Size);
};

/// @brief  Maker class that produces a new FragmentCrmsdResDepth object
class MakeFragmentCrmsdResDepth: public MakeFragmentScoringMethod {
public:

	MakeFragmentCrmsdResDepth() :
		MakeFragmentScoringMethod("FragmentCrmsdResDepth") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_FragmentCrmsdResDepth_HH
