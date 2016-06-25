// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/FragmentScoringMethod.hh
/// @brief  a base class for fragment scoring
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_FragmentScoringMethod_hh
#define INCLUDED_protocols_frag_picker_scores_FragmentScoringMethod_hh

// type headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.fwd.hh>
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <string>

namespace protocols {
namespace frag_picker {
namespace scores {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

/// @brief  a fragment candidate score
class FragmentScoringMethod: public utility::pointer::ReferenceCount {
public:

	FragmentScoringMethod(Size priority, Real lowest_acceptable_value,
		bool use_lowest, std::string name):
		priority_ (priority),
		name_(name),
		lowest_acceptable_value_(lowest_acceptable_value),
		use_lowest_(use_lowest)
	{}

	virtual bool score(FragmentCandidateOP, FragmentScoreMapOP) =0;

	/// @brief Returns a name of this scoring method
	inline std::string& get_score_name() {
		return name_;
	}

	/// @brief Returns an integer index assigned to this scoring method by a scoring manager
	inline Size get_id() {
		return id_;
	}

	/// @brief Sets a new integer index for this scoring method.
	/// @details the method should be called only by a scoring manager
	inline void set_id(Size id) {
		id_ = id;
	}

	/// @brief Returns a priority of this scoring method.
	/// @details The higher the priority, the earlier a given scoring method is evaluated
	inline Size get_priority() {
		return priority_;
	}

	/// @brief Returns the lowest acceptable score value for this scoring method.
	/// @details Fragments that are below the threshold will be discarded
	inline Real get_min_allowed_score() {
		return lowest_acceptable_value_;
	}

	/// @brief Returns the boolean choice on using the above lowest acceptable score value.
	/// @details False means there is no lowest acceptable score
	inline bool get_use_lowest() {
		return use_lowest_;
	}

	/// @brief Sets a new value of the lowest acceptable score
	/// @details Fragments that are below the threshold will be discarded
	inline void set_min_allowed_score(Real lowest_acceptable_value) {
		lowest_acceptable_value_ = lowest_acceptable_value;
	}

protected:
	Size id_;
	Size priority_;
	std::string name_;
	Real lowest_acceptable_value_;
	bool use_lowest_;
};

/// @brief  a fragment candidate
class MakeFragmentScoringMethod: public utility::pointer::ReferenceCount {
public:
	MakeFragmentScoringMethod(std::string name) :
		name_(name) {
	}

	virtual FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP,
		std::string) =0;

	inline std::string& get_score_name() {
		return name_;
	}

private:
	std::string name_;
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_FragmentScoringMethod_HH */
