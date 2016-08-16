// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/phosphate/PhosphateMove.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/phosphate/PhosphateMove.hh>
#include <map>
#include <iostream>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.phosphate.PhosphateMove" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace phosphate {

//Destructor
PhosphateMove::~PhosphateMove()
{}

///////////////////////////////////////////////
std::string
to_string( PhosphateTerminus const & phosphate_terminus ){

	static bool init( false );
	static std::map< PhosphateTerminus, std::string> phosphate_terminus_name;

	if ( !init ) {
		phosphate_terminus_name[ NONE ] = "NONE";
		phosphate_terminus_name[ FIVE_PRIME_PHOSPHATE ]  = "FIVE_PRIME_PHOSPHATE";
		phosphate_terminus_name[ THREE_PRIME_PHOSPHATE ] = "THREE_PRIME_PHOSPHATE";
		init = true;
	}

	return phosphate_terminus_name[ phosphate_terminus ];
}

/////////////////////////////////////////////////////////////////////////////////////////
bool
PhosphateMove::operator== ( PhosphateMove const & other) const{
	return ( rsd_ == other.rsd() && terminus_ == other.terminus() );
}

/////////////////////////////////////////////////////////////////////////////////////////
std::ostream &
operator <<( std::ostream & os, PhosphateMove const & phosphate_move )
{
	os << " res " << phosphate_move.rsd() << " at terminus: " << to_string( phosphate_move.terminus() );
	return os;
}

} //phosphate
} //rna
} //modeler
} //stepwise
} //protocols
