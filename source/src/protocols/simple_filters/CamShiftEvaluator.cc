// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @details
///
///
/// @author Nikolas Sgourakis ?


// Unit Headers
#include <protocols/simple_filters/CamShiftEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>


// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>

// C++ headers

static thread_local basic::Tracer tr( "protocols.simple_filter.CamShiftEvaluator" );

namespace protocols {
namespace simple_filters {

using namespace core;

#ifdef  __native_client__
#define system(a) 1
#endif


CamShiftEvaluator::CamShiftEvaluator( std::string tag, std::string cs_file )
: ExternalEvaluator( tag )
{
#ifdef WIN32
	utility_exit_with_message("don't use CamShiftEvaluator on a BillBox");
#endif
	// if (!utility::file::file_exists( scratch_dir()+"/SPARTA" ) ) {
	// std::string command( "cp -Rf $HOME/SPARTA "+scratch_dir());
	std::string command( "rsync -azvu $HOME/camshift-1.35.0 "+scratch_dir());
	int ret(system(command.c_str()));
	if ( ret ) {
		utility_exit_with_message("System command failed:'" + command + "'" );
	}


	std::string command2( "rsync -azvu $HOME/scripts/calculate_cs_rms.pl "+scratch_dir());
	int ret2(system(command2.c_str()));
	if ( ret2 ) {
		utility_exit_with_message("System command failed:'" + command2 + "'" );
	}

	//}
	set_command( "export CAMSHIFT_DIR="+scratch_dir()+
		"/camshift-1.35.0; $CAMSHIFT_DIR/bin/camshift --data ~/camshift-1.35.0/data/ --pdb __POSE.pdb >tmpc_RESULT "
		"; perl " +scratch_dir()+"/calculate_cs_rms.pl "+cs_file+" tmpc_RESULT > __RESULT" );
}

bool CamShiftEvaluator::applicable( pose::Pose const& pose ) const {
	return pose.is_fullatom();
}


}
}
