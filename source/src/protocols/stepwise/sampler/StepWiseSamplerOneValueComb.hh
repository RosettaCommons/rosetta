// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerOneValueComb.hh
/// @brief Aggregate of multiple rotamer samplers for modeler combinatorially.
/// @author  Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_protocols_sampler_StepWiseSamplerOneValueComb_HH
#define INCLUDED_protocols_sampler_StepWiseSamplerOneValueComb_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneValueComb.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerOneValueComb : public StepWiseSamplerSizedComb {

public:

	StepWiseSamplerOneValueComb();

	virtual ~StepWiseSamplerOneValueComb();

	ValueList const &
	get_value_list( Size const id );

	ValueList const &
	get_value_list( utility::vector1< Size > const & id_list );

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_external_loop_rotamer( StepWiseSamplerOneValueOP const & rotamer ) {
		StepWiseSamplerSizedComb::add_external_loop_rotamer( rotamer );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "StepWiseSamplerOneValueComb"; }

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return ONE_VALUE_COMB; }

private:

	/// @brief Add one more rotamer sampler to this sampler -- note that this is now made 'private', so that users can only add StepWiseSamplerOneValueOPs
	using	StepWiseSamplerSizedComb::add_external_loop_rotamer;

private:

	ValueList value_list_cached_, id_list_cached_;

};

} //sampler
} //stepwise
} //protocols

#endif

