// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerSized.hh
/// @brief Abstract Base Class for StepWiseSampler sampler with finite size.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_sampler_StepWiseSamplerSized_HH
#define INCLUDED_protocols_sampler_StepWiseSamplerSized_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerSized.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerSized: public StepWiseSamplerBase {

public:
	StepWiseSamplerSized();

	virtual ~StepWiseSamplerSized();

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const = 0;

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer.
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Check if reach the end of rotamer list
	virtual bool not_end() const;

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose&, core::Size const ) = 0;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose ) { apply( pose, id_ ); }

	/// @brief Name of the class
	virtual std::string get_name() const { return "StepWiseSamplerSized"; }

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return SIZED; }

	Size const & id() const { return id_; }

	virtual void fast_forward(){ id_ = size(); }

	virtual void set_id( Size const setting ){ id_ = setting; }

	/// @brief output summary of class
	virtual
	void show( std::ostream & out, Size const indent = 0 ) const {
		for ( Size n = 1; n <= indent; n++ ) out << ' ';
		out << get_name() << " [" << size() << ']' << std::endl;
	}

protected:

	core::Size id_;

};

} //sampler
} //stepwise
} //protocols
#endif
