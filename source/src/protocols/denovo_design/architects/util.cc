// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/architects/util.cc
/// @brief utility functions for architects
/// @author Tom Linsky (tlinsky@uw.edu)

#include <protocols/denovo_design/architects/util.hh>

// Protocol headers
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.util" );

namespace protocols {
namespace denovo_design {
namespace architects {

DeNovoArchitectCOP
retrieve_denovo_architect( std::string const & architect_name, basic::datacache::DataMap & data )
{
	if ( !data.has( DeNovoArchitect::DATA_MAP_NAME, architect_name ) ) {
		std::stringstream msg;
		msg << "protocols::denovo_design::architects::retrieve_denovo_architect(): architect named "
			<< architect_name << " was not found in the DataMap!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return data.get_ptr< DeNovoArchitect const >( DeNovoArchitect::DATA_MAP_NAME, architect_name );
}

void
store_denovo_architect( DeNovoArchitectOP architect, basic::datacache::DataMap & data )
{
	debug_assert( architect );
	if ( data.has( DeNovoArchitect::DATA_MAP_NAME, architect->id() ) ) {
		std::stringstream msg;
		msg << "protocols::denovo_design::architects::store_denovo_architect(): architect named "
			<< architect->id() << " already exists in the DataMap!  Architects must have unique names." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	data.add( DeNovoArchitect::DATA_MAP_NAME, architect->id(), architect );
}

} //protocols
} //denovo_design
} //architects
