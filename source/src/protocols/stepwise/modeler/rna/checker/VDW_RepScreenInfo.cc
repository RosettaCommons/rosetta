// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/checker/VDW_RepScreenInfo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/checker/VDW_RepScreenInfo.hh>
#include <core/pose/Pose.hh>
#include <string>


#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.modeler.rna.checker.VDW_RepScreenInfo" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {

//Constructor
VDW_RepScreenInfo::VDW_RepScreenInfo():
	input_string( "" ),
	pose_name( "" ),
	in_root_partition( false ),
	import_ID( 0 )
{
	VDW_align_res.clear();
	working_align_res.clear();
	full_align_res.clear();
	VDW_ignore_res.clear();
}

// Copy Constructor
VDW_RepScreenInfo::VDW_RepScreenInfo( VDW_RepScreenInfo const & src ):
	VDW_align_res( src.VDW_align_res ),
	working_align_res( src.working_align_res ),
	full_align_res( src.full_align_res ),
	VDW_ignore_res( src.VDW_ignore_res ),
	VDW_pose( src.VDW_pose ),
	input_string( src.input_string ),
	pose_name( src.pose_name ),
	in_root_partition( src.in_root_partition ),
	import_ID( src.import_ID )
{

}




//Destructor
VDW_RepScreenInfo::~VDW_RepScreenInfo()
{}

} //checker
} //rna
} //modeler
} //stepwise
} //protocols
