// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerSizedAny.hh
/// @brief Aggregate multiple samplers for modeler from any one of them.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_sampler_StepWiseSamplerSizedAny_HH
#define INCLUDED_protocols_sampler_StepWiseSamplerSizedAny_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerSizedAny.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerSizedAny : public StepWiseSamplerSized {
public:
	StepWiseSamplerSizedAny();

	virtual ~StepWiseSamplerSizedAny();

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose );

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, Size const i );

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const {
		runtime_assert( is_init() );
		return size_;
	}

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_external_loop_rotamer( StepWiseSamplerSizedOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		size_list_.clear();
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "StepWiseSamplerSizedAny"; }

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return SIZED_ANY; }

	/// @brief output summary of class
	virtual
	void show( std::ostream & out, Size const indent = 0 ) const;

private:
	/// @brief Convert an id number to the sampler state pair.
	std::pair<core::Size, core::Size> id2state( core::Size const id ) const;

	core::Size size_;

	std::pair<core::Size, core::Size> curr_state_;

	utility::vector1<core::Size> size_list_;

	utility::vector1<StepWiseSamplerSizedOP> rotamer_list_;
};

} //sampler
} //stepwise
} //protocols

#endif

