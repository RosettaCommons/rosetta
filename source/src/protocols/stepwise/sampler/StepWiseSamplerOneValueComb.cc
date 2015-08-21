// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerOneValueComb.cc
/// @brief Aggregate of multiple rotamer samplers for modeler combinatorially.
/// @author Rhiju Das (rhiju@stanford.edu)

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneValueComb.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {

///////////////////////////////////////////////////////////////////////////
StepWiseSamplerOneValueComb::StepWiseSamplerOneValueComb():
	StepWiseSamplerSizedComb()
{}

StepWiseSamplerOneValueComb::~StepWiseSamplerOneValueComb(){}


///////////////////////////////////////////////////////////////////////////
ValueList const &
StepWiseSamplerOneValueComb::get_value_list( utility::vector1< Size > const & id_list ){

	if ( id_list_cached_.size() != id_list.size() ) {
		id_list_cached_.clear();
		value_list_cached_.clear();
		for ( Size n = 1; n <= id_list.size(); n++ ) {
			id_list_cached_.push_back( 0 );
			value_list_cached_.push_back( 0.0 );
		}
	}
	runtime_assert( id_list.size() == rotamer_list_.size() );
	runtime_assert( id_list.size() == id_list_cached_.size() );

	for ( Size n = 1; n <= id_list.size(); n++ ) {
		if ( id_list[n] == id_list_cached_[n] ) continue;
		StepWiseSamplerOneValue * rotamer_one_value = static_cast< StepWiseSamplerOneValue * >( rotamer_list_[n].get() );
		value_list_cached_[n] = rotamer_one_value->value( id_list[n] );
	}

	return value_list_cached_;
}

///////////////////////////////////////////////////////////////////////////
ValueList const &
StepWiseSamplerOneValueComb::get_value_list( Size const id ){
	utility::vector1 <Size> id_list = id2list( id );
	return get_value_list( id_list );
}


} //sampler
} //stepwise
} //protocols
