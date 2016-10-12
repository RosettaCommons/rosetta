// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/AtomPairConstraintsScore.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_AtomPairConstraintsScore_hh
#define INCLUDED_protocols_frag_picker_scores_AtomPairConstraintsScore_hh

#include <protocols/frag_picker/scores/AtomPairConstraintsScore.fwd.hh>
#include <protocols/frag_picker/scores/AtomBasedConstraintsScore.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
// mini
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <numeric/xyzVector.hh>

#include <core/scoring/func/Func.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief Scores a fragment with a set of AtomPair constraints
class AtomPairConstraintsScore: public AtomBasedConstraintsScore {
public:

	AtomPairConstraintsScore(core::Size, core::Real, bool, std::string, core::Size, utility::vector1<
		std::string>);

	AtomPairConstraintsScore(core::Size, core::Real, bool, std::string, core::Size);

	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	utility::vector1<utility::vector1<AtomPairConstraintsDataOP> > data_;
	core::scoring::func::FuncFactory factory_;
	void read_constraints(std::string);
	core::Size get_atom_type(std::string atom_name);
};

/// @brief Holds data about a single distance constraint in the form AtomPairConstraintsScore needs
/// @details This class is used by AtomPairConstraintsScore to store data obtained from file
class AtomPairConstraintsData: public utility::pointer::ReferenceCount {
public:

	/// @brief makes a new object
	AtomPairConstraintsData(core::scoring::func::FuncOP function,
		core::Size offset, core::Size first_atom, core::Size second_atom) {
		offset_ = offset;
		func_ = function;
		first_atom_ = first_atom;
		second_atom_ = second_atom;
	}

	virtual ~AtomPairConstraintsData() ; // auto-removing definition from header{}

	inline core::Size get_offset() {
		return offset_;
	}

	inline core::Size get_first_atom() {
		return first_atom_;
	}

	inline core::Size get_second_atom() {
		return second_atom_;
	}

	inline core::scoring::func::FuncOP get_function() {
		return func_;
	}

private:
	core::Size offset_;
	core::scoring::func::FuncOP func_;
	core::Size first_atom_;
	core::Size second_atom_;
};

/// @brief  Maker class that produces a new AtomPairConstraintsScore object
class MakeAtomPairConstraintsScore: public MakeFragmentScoringMethod {
public:

	MakeAtomPairConstraintsScore() :
		MakeFragmentScoringMethod("AtomPairConstraintsScore") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string);
};
} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_AtomPairConstraintsScore_HH */
