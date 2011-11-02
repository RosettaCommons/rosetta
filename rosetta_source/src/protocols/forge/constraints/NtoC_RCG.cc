// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/NtoC_RCG.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009

// Unit header
#include <protocols/forge/constraints/NtoC_RCG.hh>

// Package headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/ScalarWeightedFunc.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/id/SequenceMapping.hh>

// Project headers
#include <basic/Tracer.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


static basic::Tracer TR( "protocols.forge.constraints.NtoC_RCG" );

namespace protocols{
namespace forge{
namespace constraints{

/// @brief
NtoC_RCG::NtoC_RCG():
	RemodelConstraintGenerator(),
	dist_( 11.0 ),
	coef_( 1.0 )
{}

/// @brief
NtoC_RCG::NtoC_RCG( Real const dist, Real const coef ):
	RemodelConstraintGenerator(),
	dist_( dist ),
	coef_( coef )
{}

/// @brief
NtoC_RCG::~NtoC_RCG() {}

/// @brief set weight
void
NtoC_RCG::set_weight( Real const coef )
{
	coef_ = coef;
}

/// @brief set distance of constraint
void
NtoC_RCG::set_distance( Real const dist )
{
	dist_ = dist;
}


/// @brief
void
NtoC_RCG::generate_remodel_constraints( Pose const & pose )
{
	using namespace core::scoring::constraints;

	std::string tag( "constraint_between_N_&_C_terminal_Calpha" );
	Real lb( 0.0 );
	Real ub( dist_ );
	Real sd( 1.0 );
	ScalarWeightedFuncOP cstfunc = new ScalarWeightedFunc( coef_, new BoundFunc( lb, ub, sd, tag ) );

  Size nres( pose.total_residue() );
	core::id::AtomID atom1( pose.residue_type( 1 ).atom_index( "CA" ), 1 );
	core::id::AtomID atom2( pose.residue_type( nres ).atom_index( "CA" ), nres );
	ConstraintOP const cst = new AtomPairConstraint( atom1, atom2, cstfunc );

	TR << "Constraints between N- and C- terminal: 1-" << nres << ", dist=" << dist_ << std::endl;

	this->add_constraint( cst );
} //generate constraints


} //namespace constraints
} //namespace forge
} //namespace protocols
