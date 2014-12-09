// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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

// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>
// AUTO-REMOVED #include <protocols/docking/stateless/SaneDockingProtocol.hh>
// AUTO-REMOVED #include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMoverFactory.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover.hh>
// AUTO-REMOVED #include <protocols/loops/Loop.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.hh>

// AUTO-REMOVED #include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragSet.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>


// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>

// AUTO-REMOVED #include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/func/LinearPenaltyFunction.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>


#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/ScoreMover.hh>
// AUTO-REMOVED #include <protocols/moves/CompositionMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/ConstraintSetMover.hh>

//#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

// AUTO-REMOVED #include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes


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
				using core::pose::setPoseExtraScore;
				setPoseExtraScore( pose, score_prefix_, rebuild_dist );
				break;
			}
		}
} // apply


}//namespace domain_assembly
}//namespace protocols


