// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/MC_Loop.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/sampler/MC_Loop.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.sampler.MC_Loop" );

namespace protocols {
namespace recces {
namespace sampler {

//Constructor
MC_Loop::MC_Loop():
	MC_Any()
{
	set_name( "MC_Loop" );
}

//Destructor
MC_Loop::~MC_Loop() = default;

///////////////////////////////////////////////////////////////////////////
void MC_Loop::operator++() {
	++curr_id_;
	if ( curr_id_ > rotamer_list_.size() ) curr_id_ = 1;
	++( *rotamer_list_[ curr_id_ ] );
}

///////////////////////////////////////////////////////////////////////////
void MC_Loop::show( std::ostream & out, Size const indent ) const {
	SamplerPlusPlus::show( out, indent );
	for ( Size k = 1; k <= rotamer_list_.size(); k++ ) {
		out << "Cycle " << k << " in loop: " << std::endl;
		rotamer_list_[k]->show( out, indent + 1 );
	}
}

} //sampler
} //recces
} //protocols
