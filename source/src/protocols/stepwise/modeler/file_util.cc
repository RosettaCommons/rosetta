// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/file_util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/file_util.hh>
#include <core/pose/Pose.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.file_util" );

namespace protocols {
namespace stepwise {
namespace modeler {

	///////////////////////////////////////////////////////////////////////
	std::string
	get_file_name( std::string const & silent_file, std::string const & tag )
	{
		int pos( silent_file.find( ".out" ) );
		runtime_assert( pos > -1 );
		std::string silent_file_sample( silent_file );
		silent_file_sample.replace( pos, 4, tag+".out" );
		return silent_file_sample;

	}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
remove_silent_file_if_it_exists( std::string const & silent_file){
	if ( utility::file::file_exists( silent_file ) ) {
		TR << TR.Red << "WARNING: silent_file " << silent_file << " already exists! removing..." << TR.Reset << std::endl;
		runtime_assert( std::remove( silent_file.c_str() ) == 0 );
	}
}

} //modeler
} //stepwise
} //protocols
