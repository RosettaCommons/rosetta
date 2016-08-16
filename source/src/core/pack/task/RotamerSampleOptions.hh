// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/RotamerSampleOptions.hh
/// @brief  class to describe rotamer construction header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_RotamerSampleOptions_hh
#define INCLUDED_core_pack_task_RotamerSampleOptions_hh

#include <string>

namespace core {
namespace pack {
namespace task {

/// @brief levels of extra rotamer sampling, i.e., step size and number of steps.
enum ExtraRotSample {
	NO_EXTRA_CHI_SAMPLES = 0,      //0
	EX_ONE_STDDEV,                 //1
	EX_ONE_HALF_STEP_STDDEV,       //2
	EX_TWO_FULL_STEP_STDDEVS,      //3
	EX_TWO_HALF_STEP_STDDEVS,      //4
	EX_FOUR_HALF_STEP_STDDEVS,     //5
	EX_THREE_THIRD_STEP_STDDEVS,   //6
	EX_SIX_QUARTER_STEP_STDDEVS,   //7
	ExtraRotSampleCardinality
};

ExtraRotSample
rot_sample_from_name( std::string const & name );

bool
is_rot_sample_name( std::string const & name );

const int EXTRACHI_CUTOFF_LIMIT = 18;

} // namespace task
} // namespace pack
} // namespace core

#endif
