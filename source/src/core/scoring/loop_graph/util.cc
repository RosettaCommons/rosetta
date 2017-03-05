// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/util.hh>
#include <core/scoring/loop_graph/LoopCycle.hh>
#include <core/scoring/loop_graph/evaluator/GaussianChainFuncPotentialEvaluator.hh>
#include <core/scoring/loop_graph/evaluator/SixDTransRotPotentialEvaluator.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "core.scoring.loop_graph.util" );

using namespace core::scoring::loop_graph::evaluator;

namespace core {
namespace scoring {
namespace loop_graph {

//////////////////////////////////////////////////////////////////////////////////////////////
void
get_loop_atom( Size const & res,
	core::pose::Pose const & pose,
	id::AtomID & atom_id,
	Vector & xyz,
	bool const takeoff /* as opposed to landing */
)
{
	using namespace pose::full_model_info;
	core::conformation::Residue const & rsd = get_residue( res, pose );
	std::string atom_name;
	if ( rsd.is_NA() ) {
		atom_name = takeoff ? " O3'" : " C5'";
	} else {
		runtime_assert( rsd.is_protein() );
		atom_name = takeoff ? " C  " : " N  ";
	}
	get_loop_atom( res, pose, atom_name, atom_id, xyz );
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
get_loop_atom( core::Size const & res,
	core::pose::Pose const & pose,
	std::string const & atom_name,
	id::AtomID & atom_id,
	Vector & xyz )
{
	using namespace pose::full_model_info;
	core::conformation::Residue const & rsd = get_residue( res, pose );
	atom_id = id::AtomID( rsd.atom_index( atom_name ), rsd.seqpos() );
	xyz = rsd.xyz( atom_name );
}


//////////////////////////////////////////////////////////////////////////////////////////////
LoopClosePotentialEvaluatorCOP
get_loop_close_potential( pose::Pose const & pose,
	LoopCycle const & loop_cycle,
	core::Real const & loop_fixed_cost,
	bool const use_6D_potential /* = false */ )
{

	if ( use_6D_potential &&
			loop_cycle.size() == 1  /* later, we could actually generalize to handle more complex loop cycles, if we can figure out SE(3) 6D convolutions */ ) {
		LoopClosePotentialEvaluatorCOP potential_evaluator = get_6D_trans_rot_potential_evaluator( loop_cycle, loop_fixed_cost, pose );
		if ( potential_evaluator != 0 ) return potential_evaluator;
	}

	// default mode -- derive a GaussianChainFunc that just depends on the distance between
	//  takeoff atom and landing atom. No orientation dependence.
	return GaussianChainFuncPotentialEvaluatorOP( new GaussianChainFuncPotentialEvaluator( loop_cycle, loop_fixed_cost, pose ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
LoopClosePotentialEvaluatorCOP
get_6D_trans_rot_potential_evaluator( LoopCycle const & loop_cycle,
	core::Real const & loop_fixed_cost,
	pose::Pose const & pose )
{
	using namespace pose::full_model_info;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size const & takeoff_pos( loop_cycle.loop(1).takeoff_pos() );
	Size const & landing_pos( loop_cycle.loop(1).landing_pos() );

	// right now, set up for RNA. Later can generalize to protein, DNA, etc.
	std::string tag;
	core::conformation::Residue const & rsd = get_residue( takeoff_pos, pose );
	if ( rsd.is_RNA() ) {
		runtime_assert( get_residue( landing_pos, pose ).is_RNA() );
		tag = "rna";
	}

	// later could generalize for cyclization (in which case takeoff_pos might actually be *after* landing_pos.
	runtime_assert( landing_pos - takeoff_pos >= 1 );
	Size loop_length( landing_pos - takeoff_pos - 1 );

	std::string database_file = basic::database::full_name( "scoring/loop_close/6D_potentials/"+tag+"/loop_"+ObjexxFCL::lead_zero_string_of( loop_length, 2 )+"/potential.txt.gz" );
	if ( option[ score::loop_close::force_6D_potential_file ].user() ) {
		database_file = option[ score::loop_close::force_6D_potential_file ]();
		runtime_assert( ScoringManager::get_instance()->get_LoopCloseSixDPotential( database_file ) != 0 );
	}
	if ( ScoringManager::get_instance()->get_LoopCloseSixDPotential( database_file ) == 0 ) return 0;

	return SixDTransRotPotentialEvaluatorCOP( new SixDTransRotPotentialEvaluator(
		takeoff_pos, landing_pos,
		pose /* needed to check if takeoff/landing of loop is in current pose */ ,
		loop_fixed_cost,
		*( ScoringManager::get_instance()->get_LoopCloseSixDPotential( database_file ) ) ) );
}





} //loop_graph
} //scoring
} //core
