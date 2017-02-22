// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/GaussianChainFuncPotentialEvaluator.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/evaluator/GaussianChainFuncPotentialEvaluator.hh>
#include <core/scoring/loop_graph/LoopCycle.hh>
#include <core/scoring/loop_graph/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/func/GaussianChainFunc.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.loop_graph.evaluator.GaussianChainFuncPotentialEvaluator" );

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

	//Constructor
	GaussianChainFuncPotentialEvaluator::GaussianChainFuncPotentialEvaluator( LoopCycle const & loop_cycle,
																																						core::Real const & loop_fixed_cost,
																																						pose::Pose const & pose ):
		LoopClosePotentialEvaluator( loop_fixed_cost ),
		rna_gaussian_variance_per_residue_( 5.0 * 5.0 ), // in Angstroms^2
		protein_gaussian_variance_per_residue_( 3.0 * 3.0 ) // in Angstroms^2
	{
		initialize( loop_cycle, pose );
	}

	//Destructor
	GaussianChainFuncPotentialEvaluator::~GaussianChainFuncPotentialEvaluator()
	{}

	void
	GaussianChainFuncPotentialEvaluator::initialize( LoopCycle const & loop_cycle,
																									 pose::Pose const & pose )
	{
		using namespace core::scoring::func;
		using namespace core::pose::full_model_info;
		using namespace core::id;
		// check list of loops, and add up total gaussian variance
		Real total_gaussian_variance( 0.0 );
		for ( Size n = 1; n <= loop_cycle.size(); n++ ) {

			Loop const & loop = loop_cycle.loop( n );
			Size const loop_length = ( loop.landing_pos() - loop.takeoff_pos() );

			if ( get_residue( loop.takeoff_pos(), pose ).is_NA() ) {
				runtime_assert( get_residue( loop.landing_pos(), pose ).is_NA() );
				total_gaussian_variance += Real( loop_length ) * rna_gaussian_variance_per_residue_;
			} else if ( get_residue( loop.takeoff_pos(), pose ).is_protein() ) {
				runtime_assert( get_residue( loop.landing_pos(), pose ).is_protein() );
				total_gaussian_variance += Real( loop_length ) * protein_gaussian_variance_per_residue_;
			} else {
				utility_exit_with_message( "loop_close GaussianChainFuncPotentialEvaluator cannot currently assign energies to loop-cycles with non-protein or non-nucleic acid parts" );
			}
		}

		// go into each 'domain' and get list of distances. Also figure out if any of the domains are the current pose, and keep track of AtomID info.
		Size current_pose_idx_in_cycle( 0 );
		AtomID takeoff_atom_id, landing_atom_id, current_pose_takeoff_atom_id, current_pose_landing_atom_id;
		Vector takeoff_xyz, landing_xyz;
		utility::vector1< Real > all_distances, other_distances;

		for ( Size n = 1; n <= loop_cycle.size(); n++ ) {
			Loop const & takeoff_loop = loop_cycle.loop( n );
			Size const & takeoff_pos  = takeoff_loop.takeoff_pos();
			Size const & takeoff_domain = takeoff_loop.takeoff_domain();

			get_loop_atom( takeoff_pos, pose, takeoff_atom_id, takeoff_xyz, true /*takeoff*/ );

			// Now need to follow cycle around to find loop that ends on this same domain. There should be only one.
			Size const landing_loop_idx  = loop_cycle.find_index_for_loop_landing_at_domain( takeoff_domain );
			Loop const & landing_loop = loop_cycle.loop( landing_loop_idx );
			Size const & landing_pos  = landing_loop.landing_pos();
			Size const & landing_domain = landing_loop.landing_domain();
			runtime_assert( takeoff_domain == landing_domain );

			get_loop_atom( landing_pos, pose, landing_atom_id, landing_xyz, false /*takeoff*/ );

			Distance d = ( landing_xyz - takeoff_xyz ).length();
			all_distances.push_back( d );
			if ( takeoff_domain == 1 ) {
				current_pose_takeoff_atom_id = takeoff_atom_id;
				current_pose_landing_atom_id = landing_atom_id;
				current_pose_idx_in_cycle = n;
			}
		}

		Size main_pose_idx = 1;
		if ( current_pose_idx_in_cycle > 0 ) main_pose_idx = current_pose_idx_in_cycle;

		// reorder so that current_distance is in the beginning.
		Real main_distance = all_distances[ main_pose_idx ];
		for ( Size n = 1; n <= all_distances.size(); n++ ) {
			if ( n != main_pose_idx ) other_distances.push_back( all_distances[ n ] );
		}

		// Note loop_fixed_cost_ needs to be my current loop_fixed_cost_ but corrected by 1.5 * k_B_T_ * log( rna_persistence_length2_ )
		func_ = FuncOP( new GaussianChainFunc( total_gaussian_variance, loop_fixed_cost(), other_distances ) );
		set_loop_closure_energy( func_->func( main_distance ) );

		if ( current_pose_idx_in_cycle ) {
			set_involves_current_pose( true );
			set_current_pose_takeoff_atom( current_pose_takeoff_atom_id );
			set_current_pose_landing_atom( current_pose_landing_atom_id );
			set_current_pose_takeoff_xyz( pose.xyz( current_pose_takeoff_atom_id ) );
			set_current_pose_landing_xyz( pose.xyz( current_pose_landing_atom_id ) );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	GaussianChainFuncPotentialEvaluator::get_f1_f2( Vector & f1, Vector & f2, bool const takeoff ) const
	{
		if ( takeoff ) {
			Vector const d = current_pose_takeoff_xyz() - current_pose_landing_xyz();
			Real const x = d.length();
			f2 = func_->dfunc( x ) * d / x;
			f1 = cross( f2,  current_pose_takeoff_xyz() );
		} else {
			// actually could just calculate f1 & f2 above, and supply negative for landing & takeoff.
			// should get same answer, up to very tiny numerical errors.
			// to maintain integration tests, let's use historical order:
			Vector const d = current_pose_landing_xyz() - current_pose_takeoff_xyz();
			Real const x = d.length();
			f2 = func_->dfunc( x ) * d / x;
			f1 = cross( f2,  current_pose_landing_xyz() );
		}
	}


} //evaluator
} //loop_graph
} //scoring
} //core
