/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/AddAssemblyConstraints.cc
/// @author James Thompson

#include <protocols/domain_assembly/AddAssemblyConstraints.hh>

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

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
#include <core/pose/PDB_Info.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>


// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>

#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>


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

#include <core/kinematics/Jump.hh>


namespace protocols {
namespace domain_assembly {

void AddAssemblyConstraints::apply( core::pose::Pose & pose ) {

		using core::Real;
		using std::string;
		using core::id::AtomID;
		string const atom_name("CA");
		using namespace core::scoring::constraints;
		Real const dist_cutoff(3), range(1.5), slope(3.0),
			score(-40);
			//score( -1 * pose.total_residue() * 0.25 );

		ConstraintSetOP cst_set( new ConstraintSet );
		char const first_chain( pose.pdb_info()->chain(1) );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			char const chain_ii( pose.pdb_info()->chain(ii) );
			if ( first_chain != chain_ii ) {
				AtomID atom1( pose.residue_type(ii-1).atom_index(atom_name), ii-1 ),
					atom2( pose.residue_type(ii).atom_index(atom_name), ii );
				core::scoring::func::FuncOP func( new core::scoring::func::LinearPenaltyFunction( dist_cutoff, score, range, slope ) );

				ConstraintOP cst = new AtomPairConstraint(atom1,atom2,func);
				cst_set->add_constraint(cst);
				break;
			}
		}
		pose.constraint_set(cst_set);
} // apply

}//namespace domain_assembly
}//namespace protocols
