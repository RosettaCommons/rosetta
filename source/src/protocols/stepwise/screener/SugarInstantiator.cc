// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/SugarInstantiator.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/SugarInstantiator.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarInstantiateMover.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarVirtualizeMover.hh>
#include <protocols/moves/CompositionMover.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.SugarInstantiator" );

using namespace protocols::stepwise::modeler::rna::sugar;

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
	SugarInstantiator::SugarInstantiator( pose::Pose & screening_pose,
																				Size const moving_res,
																				Distance const o2prime_instantiation_distance_cutoff ):
		screening_pose_( screening_pose ),
		moving_res_( moving_res ),
		o2prime_instantiation_distance_cutoff_( o2prime_instantiation_distance_cutoff ),
		instantiate_sugar_( false )
	{}

	//Destructor
	SugarInstantiator::~SugarInstantiator()
	{}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	SugarInstantiator::check_screen(){
		instantiate_sugar_ = check_moving_sugar( screening_pose_, moving_res_ );
		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	SugarInstantiator::check_moving_sugar( pose::Pose & pose, Size const moving_res ){

		if ( !pose.residue( moving_res_ ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) return false; // nothing to do.

		bool instantiate_sugar( false );
		Vector const & moving_O2prime_xyz = pose.residue( moving_res ).xyz( " O2'" );
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			if ( i == moving_res ) continue;
			if ( pose.residue( i ).is_virtual_residue() ) continue;
			for ( Size j = 1; j <= pose.residue( i ).nheavyatoms(); j++ ){
				if ( pose.residue_type( i ).is_virtual( j ) ) continue;
				if ( pose.residue_type( i ).heavyatom_is_an_acceptor( j ) ||
						 pose.residue_type( i ).heavyatom_has_polar_hydrogens( j )  ){
					Distance dist = ( pose.residue(i).xyz( j ) - moving_O2prime_xyz ).length();
					//std::cout << "Distance to rsd " << i << "  atom " << pose.residue_type( i ).atom_name( j ) << ": " << dist << std::endl;
					if ( dist < o2prime_instantiation_distance_cutoff_ ) {
						instantiate_sugar = true; break;
					}
				}
			} // atoms j
			if ( instantiate_sugar ) break;
		} // residues i

		return instantiate_sugar;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SugarInstantiator::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
		if ( instantiate_sugar_ ){
			using protocols::moves::MoverOP;
			update_mover->add_mover( new SugarInstantiateMover( moving_res_) );
			restore_mover->add_mover( new SugarVirtualizeMover( moving_res_) );
		} else {
			update_mover->add_mover( 0 );
			restore_mover->add_mover( 0 );
		}
	}



} //screener
} //stepwise
} //protocols
