// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/RotamerSampleOptions.cc
/// @brief
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#include <string>

#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/types.hh>

#include <utility/exit.hh>

#include <utility/vector0_bool.hh>


namespace core {
namespace pack {
namespace task {

ExtraRotSample
rot_sample_from_name( std::string const & name )
{
	utility::vector0< std::string > rotsample_names;
	rotsample_names.reserve( ExtraRotSampleCardinality );
	rotsample_names.push_back( "NO_EXTRA_CHI_SAMPLES" );
	rotsample_names.push_back( "EX_ONE_STDDEV" );
	rotsample_names.push_back( "EX_ONE_HALF_STEP_STDDEV" );
	rotsample_names.push_back( "EX_TWO_FULL_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_TWO_HALF_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_FOUR_HALF_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_THREE_THIRD_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_SIX_QUARTER_STEP_STDDEVS" );

	for ( Size ii = 0; ii < ExtraRotSampleCardinality; ++ii ) {
		if ( name == rotsample_names[ ii ] ) {
			return ExtraRotSample( ii );
		}
	}

	utility_exit_with_message( "ERROR: Could not find ExtraRotSample ID for string '" + name + "'" );
	return NO_EXTRA_CHI_SAMPLES;
}


bool
is_rot_sample_name( std::string const & name )
{
	utility::vector0< std::string > rotsample_names;
	rotsample_names.reserve( ExtraRotSampleCardinality );
	rotsample_names.push_back( "NO_EXTRA_CHI_SAMPLES" );
	rotsample_names.push_back( "EX_ONE_STDDEV" );
	rotsample_names.push_back( "EX_ONE_HALF_STEP_STDDEV" );
	rotsample_names.push_back( "EX_TWO_FULL_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_TWO_HALF_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_FOUR_HALF_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_THREE_THIRD_STEP_STDDEVS" );
	rotsample_names.push_back( "EX_SIX_QUARTER_STEP_STDDEVS" );

	for ( Size ii = 0; ii < ExtraRotSampleCardinality; ++ii ) {
		if ( name == rotsample_names[ ii ] ) {
			return true;
		}
	}
	return false;
}


} // namespace task
} // namespace pack
} // namespace core

