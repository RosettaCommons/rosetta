// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/DebuggingAnnealer.cc
/// @brief  Debugging annealer class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/DebuggingAnnealer.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

#include <ObjexxFCL/Fmath.hh>

#include <utility/io/izstream.hh>

#include <fstream>
#include <istream>
#include <iostream>

namespace core {
namespace pack {
namespace annealer {

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// constructor
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author


////////////////////////////////////////////////////////////////////////////////
DebuggingAnnealer::DebuggingAnnealer(
	utility::vector0< int > & rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::InteractionGraphBaseOP ig,
	rotamer_set::FixbbRotamerSetsCOP p_rotamer_set,
	ObjexxFCL::FArray1_int & current_rot_index,
	bool calc_rot_freq,
	ObjexxFCL::FArray1D_float & rot_freq
):
	RotamerAssigningAnnealer(
	rot_to_pack,
	(int) rot_to_pack.size(),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	p_rotamer_set,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	), ig_(ig)
{
}

DebuggingAnnealer::DebuggingAnnealer(
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::InteractionGraphBaseOP ig,
	rotamer_set::FixbbRotamerSetsCOP p_rotamer_set,
	ObjexxFCL::FArray1_int & current_rot_index,
	bool calc_rot_freq,
	ObjexxFCL::FArray1D_float & rot_freq
):
	RotamerAssigningAnnealer(
	(ig->get_num_total_states()),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	p_rotamer_set,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	), ig_(ig)
{
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// virtual destructor
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
DebuggingAnnealer::~DebuggingAnnealer()
{}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// sim_annealing for fixed backbone design mode
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
void DebuggingAnnealer::run()
{

	//--------------------------------------------------------------------
	//internal variables

	int const nmoltenres = ig_->get_num_nodes();
	//float currentenergy;  // unused ~Labonte
	//FArray1D_int list( p_rotamer_set_->nrotamers() );
	//FArray1D_int state_on_node( nmoltenres,0 );
	ObjexxFCL::FArray1D_int best_state_on_node( nmoltenres,0 );
	ObjexxFCL::FArray1D_int current_state( nmoltenres, 0 );
	//FArray1D_float loopenergy(maxouteriterations,0.0);

	//--------------------------------------------------------------------
	//initialize variables

	//currentenergy = 0.0;  // unused ~Labonte

	ig_->prepare_for_simulated_annealing();
	ig_->blanket_assign_state_0();

	//--------------------------------------------------------------------
	if ( num_rots_to_pack() == 0 ) return;

	setup_iterations();

	ObjexxFCL::FArray1D_float previous_nsteps_for_rot(rotamer_sets()->nrotamers(), 0.0);

	//int outeriterations = get_outeriterations();


	for ( std::list< RotSub >::const_iterator
			iter = trajectory_.begin(), iter_end = trajectory_.end();
			iter != iter_end; ++iter ) {

		int node_to_change = iter->moltenresid;
		int new_state_for_node = iter->rotamerid;

		float deltaE( 0.0f ), previous_energy_for_node( 0.0f );
		ig_->consider_substitution( node_to_change, new_state_for_node, deltaE, previous_energy_for_node );

		//std::cout << "mres: " << node_to_change << ", state: ";
		//std::cout << new_state_for_node << ", deltaE: " << deltaE << std::endl;

		if ( iter->accept ) {
			ig_->commit_considered_substitution();
			current_state( node_to_change ) = new_state_for_node;
		}

	}

	best_state_on_node = current_state;

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( int ii = 1; ii <= nmoltenres; ++ii ) {
		int iiresid = rotamer_sets()->moltenres_2_resid(ii);
		bestrotamer_at_seqpos()(iiresid) = best_state_on_node(ii) + rotamer_sets()->nrotamer_offset_for_moltenres(iiresid);
	}

}

void DebuggingAnnealer::annealer_file( std::string const & fname )
{
	trajectory_.clear();
	utility::io::izstream instructions( fname.c_str() );
	while ( instructions ) {
		RotSub rotsub;
		instructions >> rotsub.moltenresid;
		if ( rotsub.moltenresid <= 0 || rotsub.moltenresid > ig_->get_num_nodes() ) {
			std::cerr << "Error in trajectory for Debugging Annealer: node_to_change out of range." << std::endl;
			std::cerr << "Requires 0 < node_in_range <= nmoltenres.  nmoltenres: " << ig_->get_num_nodes() << " node_to_change: " << rotsub.moltenresid << std::endl;
			return;
		}

		instructions >> rotsub.rotamerid;

		if ( rotsub.rotamerid <= 0 || rotsub.rotamerid > ig_->get_num_states_for_node( rotsub.moltenresid ) ) {
			std::cerr << "Error in instructions for Debugging Annealer: new_state_for_node out of range." << std::endl;
			std::cerr << "Requires 0 < new_state_for_node <= num_states_for_node( " <<  rotsub.moltenresid << ");" << std::endl;
			std::cerr << "new_state_for_node: " << rotsub.rotamerid << " num_states_for_node( " <<  rotsub.moltenresid << "): ";
			std::cerr << ig_->get_num_states_for_node( rotsub.moltenresid ) << std::endl;
			return;
		}

		char decision;
		instructions >> decision;
		if ( decision == 'A' ) {
			rotsub.accept = 1;
		} else {
			rotsub.accept = 0;
		}
		trajectory_.push_back( rotsub );
	}
}


}//end namespace annealer
}//end namespace pack
}//end namespace core
