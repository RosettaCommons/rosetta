// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/domain_assembly/PostDockAssemblyScorer.cc
/// @brief  Computes crmsd of the assembly to the staring point
/// @author James Thompson

#include <protocols/domain_assembly/PostDockAssemblyScorer.hh>


#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>


#include <protocols/moves/Mover.hh>

//#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes


namespace protocols {
namespace domain_assembly {

void PostDockAssemblyScorer::apply( core::pose::Pose & pose ) {
	using core::Real;
	using std::string;

	char const first_chain( pose.pdb_info()->chain(1) );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		char const chain_ii( pose.pdb_info()->chain(ii) );
		if ( first_chain != chain_ii ) {
			Real const rebuild_dist(
				pose.residue(ii-1).xyz("CA").distance(
				pose.residue(ii).xyz("CA")
				)
			);
			using core::pose::setPoseExtraScore;
			setPoseExtraScore( pose, score_prefix_, rebuild_dist );
			break;
		}
	}
} // apply


}//namespace domain_assembly
}//namespace protocols


