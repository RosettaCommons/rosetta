// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/InterbondAngleScore.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_InterbondAngleScore_hh
#define INCLUDED_protocols_frag_picker_scores_InterbondAngleScore_hh

#include <protocols/frag_picker/scores/InterbondAngleScore.fwd.hh>
#include <protocols/frag_picker/scores/AtomBasedConstraintsScore.hh>

// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
// mini
#include <core/scoring/func/FuncFactory.hh>
#include <numeric/xyzVector.hh>

#include <protocols/frag_picker/scores/FourAtomsConstraintData.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief Scores a fragment with a set of Dihedral constraints
class InterbondAngleScore: public AtomBasedConstraintsScore {
public:

	/// @brief Prepare an atom-based score that utilizes some user-defined atoms
	/// @details User may provide names of atoms that will be cached when a new
	/// chunk is considered (i.e. at every do_caching() call)
	InterbondAngleScore(core::Size, core::Real, bool, std::string, core::Size, utility::vector1<
		std::string>);

	/// @brief Prepare an atom-based score that utilizes the following predefined atoms: N, CA, C, O and CB
	/// @details These atoms that will be cached when a new
	/// chunk is considered (i.e. at every do_caching() call)
	InterbondAngleScore(core::Size, core::Real, bool, std::string, core::Size);

	/// @brief Calculates the score
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	utility::vector1<utility::vector1<FourAtomsConstraintDataOP> > data_;
	void read_constraints(std::string);
	core::Size get_atom_type(std::string atom_name);
	core::scoring::func::FuncFactory factory_;
};


/// @brief  Maker class that produces a new InterbondAngleScore object
class MakeInterbondAngleScore: public MakeFragmentScoringMethod {
public:

	MakeInterbondAngleScore() :
		MakeFragmentScoringMethod("InterbondAngleScore") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string);
};
} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_InterbondAngleScore_HH */
