// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/domain_assembly/AddAssemblyConstraints.cc
/// @author James Thompson

#include <protocols/domain_assembly/AddAssemblyConstraints.hh>

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>


#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>


#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>


//#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

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
	//score( -1 * pose.size() * 0.25 );

	ConstraintSetOP cst_set( new ConstraintSet );
	char const first_chain( pose.pdb_info()->chain(1) );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		char const chain_ii( pose.pdb_info()->chain(ii) );
		if ( first_chain != chain_ii ) {
			AtomID atom1( pose.residue_type(ii-1).atom_index(atom_name), ii-1 ),
				atom2( pose.residue_type(ii).atom_index(atom_name), ii );
			core::scoring::func::FuncOP func( new core::scoring::func::LinearPenaltyFunction( dist_cutoff, score, range, slope ) );

			ConstraintOP cst( new AtomPairConstraint(atom1,atom2,func) );
			cst_set->add_constraint(cst);
			break;
		}
	}
	pose.constraint_set(cst_set);
} // apply

}//namespace domain_assembly
}//namespace protocols
