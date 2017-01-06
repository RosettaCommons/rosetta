// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/MC_Any.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


// Unit headers
#include <protocols/recces/sampler/MC_Any.hh>

// Project headers
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static basic::Tracer TR( "protocols.recces.sampler.MC_Any" );

namespace protocols {
namespace recces {
namespace sampler {

	///////////////////////////////////////////////////////////////////////////
	//Constructor
	MC_Any::MC_Any():
		MC_Comb()
	{
		set_name( "MC_Any" );
	}
	///////////////////////////////////////////////////////////////////////////
	//Destructor
	MC_Any::~MC_Any()
	{}
	///////////////////////////////////////////////////////////////////////////
	void MC_Any::init() {
		MC_Comb::init();
		curr_id_ = 0;
	}
	///////////////////////////////////////////////////////////////////////////
	void MC_Any::operator++() {
		curr_id_ = numeric::random::rg().random_range(1,num_rotamers());
		++( *rotamer_list_[ curr_id_ ] );
	}
	///////////////////////////////////////////////////////////////////////////
	void MC_Any::apply( Pose & pose ) {
		runtime_assert( is_init() );
		rotamer_list_[ curr_id_ ]->apply( pose );
		found_move_ = rotamer_list_[ curr_id_ ]->found_move();
	}

} //sampler
} //recces
} //protocols
