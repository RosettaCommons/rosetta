/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/PostDockAssemblyScorer.cc
/// @brief  Computes crmsd of the assembly to the staring point
/// @author James Thompson

#include <protocols/domain_assembly/PostDockAssemblyScorer.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/docking/stateless/SaneDockingProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/loops/loops_main.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/LinearPenaltyFunction.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>


#include <protocols/moves/Mover.hh>
#include <protocols/moves/ScoreMover.hh>
#include <protocols/moves/CompositionMover.hh>
#include <protocols/moves/ConstraintSetMover.hh>

//#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes

#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers
#include <basic/options/option.hh>


namespace protocols {
namespace domain_assembly {

void PostDockAssemblyScorer::apply( core::pose::Pose & pose ) {
		using core::Real;
		using std::string;

		char const first_chain( pose.pdb_info()->chain(1) );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			char const chain_ii( pose.pdb_info()->chain(ii) );
			if ( first_chain != chain_ii ) {
				Real const rebuild_dist(
					pose.residue(ii-1).xyz("CA").distance(
						pose.residue(ii).xyz("CA")
					)
				);
				using core::pose::setPoseExtraScores;
				setPoseExtraScores( pose, score_prefix_, rebuild_dist );
				break;
			}
		}
} // apply


}//namespace domain_assembly
}//namespace protocols


