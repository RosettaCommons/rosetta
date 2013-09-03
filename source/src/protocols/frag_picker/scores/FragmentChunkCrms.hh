// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/FragmentChunkCrms.hh
/// @brief  Object that scores a fragment by its tmscore to the native
/// @author Lei Shi (shilei@uw.edu)

#ifndef INCLUDED_protocols_frag_picker_scores_FragmentChunkCrms_hh
#define INCLUDED_protocols_frag_picker_scores_FragmentChunkCrms_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// mini
// AUTO-REMOVED #include <protocols/toolbox/superimpose.hh>

#include <core/types.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

//Auto Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <iostream>

#include <math.h>
#include <vector>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <core/sequence/Sequence.fwd.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

using namespace ObjexxFCL;

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  scores a fragment by its tmscore to the given reference structure
class FragmentChunkCrms: public FragmentScoringMethod {
public:

	/// @brief  creates a tmscore-based scoring function.
	/// @detailed fragments will be compared to a given pose, which should have the same number of residues a the query sequence
	FragmentChunkCrms(Size, Real, bool, std::string, core::pose::PoseOP, FArray1D_int&);

	~FragmentChunkCrms();

	bool score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	Size n_atoms_;
  //std::string query_sequence_;
	pose::PoseOP reference_pose_;
  pose::PoseOP fragment_pose_;
	FArray2D_double reference_coordinates_;
	FArray2D_double chunk_coordinates_;
	FArray2D_double fragment_coordinates_;
	FArray1D_double weights_;
	FArray1D_int seqmapping_;

	void fill_bb_coords(core::pose::Pose const &, FArray2_double &, Size);
	void fill_bb_coords(core::pose::Pose const &, FArray2_double &, FArray1D_int &);
};

// Undefined, commenting out to fix PyRosetta build
// void sequencealign(core::sequence::SequenceOP seq1, core::sequence::SequenceOP seq2);

/// @brief  Maker class that produces a new FragmentChunkCrms object
class MakeFragmentChunkCrms: public MakeFragmentScoringMethod {
public:

	MakeFragmentChunkCrms() :
		MakeFragmentScoringMethod("FragmentChunkCrms") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_FragmentChunkCrms_HH
