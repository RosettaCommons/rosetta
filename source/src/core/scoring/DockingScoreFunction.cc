// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  ScoreFunction class definition.
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Modified by Sergey Lyskov


// Unit headers
#include <core/scoring/DockingScoreFunction.hh>

// Package headers

// // Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer tr( "core.scoring.DockingScoreFunction" );

namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////
DockingScoreFunction::DockingScoreFunction():
	ScoreFunction()
{}

///////////////////////////////////////////////////////////////////////////////
ScoreFunctionOP
DockingScoreFunction::clone() const
{
	DockingScoreFunctionOP newscorefxn( new DockingScoreFunction );
	newscorefxn->assign( *this );
	return newscorefxn;
}

///////////////////////////////////////////////////////////////////////////////

/// @brief INTERNAL USE ONLY
void
DockingScoreFunction::assign( DockingScoreFunction const & src )
{
	ScoreFunction::assign( src );
	// Add assignemnt of DockingScoreFunction specific values here.
}

/// @brief INTERNAL USE ONLY
void
DockingScoreFunction::assign( ScoreFunction const & src )
{
	ScoreFunction::assign( src );
}

///////////////////////////////////////////////////////////////////////////////

// to start out, just thinking fullatom energies
//
// NOTE: no freakin rotamer trials inside scoring!
Real
DockingScoreFunction::operator()( pose::Pose & pose ) const
{
	Size interface_jump_id = 1; //for now

	core::kinematics::Jump bound_pose_jump = pose.jump( interface_jump_id );
	core::kinematics::Jump unbound_pose_jump( bound_pose_jump );
	Real trans_magnitude = 10000;

	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( interface_jump_id );
	core::Vector dummy_axis(1,1,1);
	unbound_pose_jump.translation_along_axis( upstream_stub, dummy_axis, trans_magnitude );
	pose.set_jump( interface_jump_id, unbound_pose_jump );

	ScoreFunction::operator()( pose ); //score -- but without atom_pair_constraints..
	// is probably cheaper to not apply a completely new scorefunction...
	EnergyMap cst_free_weights( weights() );
	Real cst_weight = cst_free_weights.get( atom_pair_constraint );
	cst_free_weights[ atom_pair_constraint ] = 0;
	Real unbound_energy = pose.energies().total_energies().dot( cst_free_weights );


	pose.set_jump( interface_jump_id, bound_pose_jump );

	ScoreFunction::operator() ( pose );
	Real bound_energy = pose.energies().total_energies().dot( cst_free_weights );

	runtime_assert( pose.num_jump() > 0 );

	Real interaction_energy = (bound_energy - unbound_energy);

	pose.energies().total_energies()[ total_score ] = interaction_energy + pose.energies().total_energies()[ atom_pair_constraint ]*cst_weight;
	tr.Debug << "unbound_energy " << unbound_energy << " full_score " << bound_energy << " diff " << interaction_energy << std::endl;
	pose::setPoseExtraScore( pose, "I_sc", interaction_energy );
	return pose.energies().total_energies()[ total_score ];
}


///////////////////////////////////////////////////////////////////////////////
} // namespace scoring
} // namespace core
