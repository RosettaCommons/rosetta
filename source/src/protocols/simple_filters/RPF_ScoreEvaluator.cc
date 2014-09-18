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
/// @detailed
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/simple_filters/RPF_ScoreEvaluator.hh>
//#include <basic/options/option.hh>

// Package Headers
#include <protocols/noesy_assign/NoesyModule.hh>
#include <protocols/noesy_assign/util.hh>
// Project Headers
#include <core/pose/Pose.hh>

// ObjexxFCL Headers

// Utility headers
//#include <basic/Tracer.hh>
//#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/util.hh>
//#include <utility/file/file_sys_util.hh>
//#include <numeric/random/random.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <iostream>
// option key includes

//#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers



// C++ headers

//static basic::Tracer tr("protocols\.simple_filter\.RPF_ScoreEvaluator");


namespace protocols {
namespace simple_filters {

using namespace core;
using namespace std;
using namespace noesy_assign;

RPF_ScoreEvaluator::RPF_ScoreEvaluator( std::string tag, core::Real dcut  )
  : evaluation::SingleValuePoseEvaluator<core::Real>( tag ),
		crosspeaks_( NULL ),
		dcut_( dcut )
{


}

bool RPF_ScoreEvaluator::applicable(  core::pose::Pose const& pose ) const {
	return pose.is_fullatom();
}

core::Real RPF_ScoreEvaluator::apply( core::pose::Pose& pose ) const {
	if ( !noesy_assign::NoesyModule::cmdline_options_activated() ) return -9999;
	if ( !crosspeaks_ ) {
		crosspeaks_ = new CrossPeakList( noesy_assign::NoesyModule( pose.sequence() ).crosspeaks() );
		crosspeaks_->find_assignments();
	}

	//return protocols::noesy_assign::compute_RPF_score( *crosspeaks_, pose, dcut_ );
}


// core::Real RPF_ScoreEvaluator::apply( core::pose::Pose& pose_in  ) const {
// 	pose::Pose pose( pose_in );

// 	runtime_assert( constraints_ );
// 	pose.constraint_set( constraints_ );

// 	ScoreFunction scfxn;
// 	scfxn.set_weight( atom_pair_constraint, 1.0 );
// 	return scfxn( pose );

// }


}
}
