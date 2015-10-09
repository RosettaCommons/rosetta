// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FixbbSimAnnealer.cc
/// @brief  Packer's standard simulated annealing class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/FixbbSimAnnealer.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <basic/Tracer.hh>

//#include "after_opts.h"
//#include "FixbbSimAnnealer.h"
//#include "RotamerAssigningAnnealer.h"
//#include "random_numbers.h"
//#include "param.h"
//#include "RotamerSet.h"

#include <utility/exit.hh>


#include <iostream>
#include <fstream>

#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace ObjexxFCL;

#ifndef NDEBUG
static THREAD_LOCAL basic::Tracer TR( "core.pack.annealer.FixbbSimAnnealer" );
#endif

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
FixbbSimAnnealer::FixbbSimAnnealer(
	utility::vector0< int > & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_sets,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D_float & rot_freq
):
	RotamerAssigningAnnealer(
	rot_to_pack,
	(int) rot_to_pack.size(),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	rotamer_sets,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	),
	ig_(ig),
	record_annealer_trajectory_( false )
{
}

FixbbSimAnnealer::FixbbSimAnnealer(
	FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D_float & rot_freq
):
	RotamerAssigningAnnealer(
	(ig->get_num_total_states()),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	rotamer_set,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	),
	ig_(ig),
	record_annealer_trajectory_( false )
{
}

/// @brief virtual destructor
FixbbSimAnnealer::~FixbbSimAnnealer()
{}

/// @brief sim_annealing for fixed backbone design mode
void FixbbSimAnnealer::run()
{

	int const nmoltenres = ig_->get_num_nodes();

	FArray1D_int state_on_node( nmoltenres,0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D_float loopenergy(maxouteriterations,0.0);

	//bk variables for calculating rotamer frequencies during simulation
	int nsteps = 0;
	FArray1D_int nsteps_for_rot( ig_->get_num_total_states(), 0 );

	//--------------------------------------------------------------------
	//initialize variables

	core::PackerEnergy currentenergy = 0.0;

	ig_->prepare_for_simulated_annealing();
	ig_->blanket_assign_state_0();

	//--------------------------------------------------------------------
	if ( num_rots_to_pack() == 0 ) return;

	setup_iterations();

	FArray1D_float previous_nsteps_for_rot( rotamer_sets()->nrotamers(), 0.0);

	int outeriterations = get_outeriterations();


	std::ofstream annealer_trajectory;
	if ( record_annealer_trajectory_ ) {
		annealer_trajectory.open(trajectory_file_name_.c_str() );
	}

	//std::cout << "Annealing begins" << std::endl;

	//outer loop
	for ( int nn = 1; nn <= outeriterations; ++nn ) {
		setup_temperature(loopenergy,nn);
		if ( quench() ) {
			currentenergy = bestenergy();
			state_on_node = best_state_on_node;
			ig_->set_network_state( state_on_node );
		}
		//rh std::cout << "Sim Annealer Temperature: " << get_temperature() << std::endl;

		int inneriterations = get_inneriterations();

		float treshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( treshold_for_deltaE_inaccuracy );

		//inner loop
		for ( int n = 1; n <= inneriterations; ++n ) {
			int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int const rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// initializing to zero but should be updated below.
			core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );

			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
				delta_energy, previous_energy_for_node);
			//std::cout << "mres: " << moltenres_id << ", state: ";
			//std::cout << rotamer_state_on_moltenres << ", deltaE: " << delta_energy;

			//bk keep new rotamer if it is lower in energy or accept it at some
			//bk probability if it is higher in energy, if it is the first
			//bk rotamer to be tried at this position automatically accept it.
			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
				//std::cout << " accepted\n";
				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					best_state_on_node = state_on_node;

					/*
					//ronj - debugging output useful for seeing how the energy changes during the course of the simulation
					for ( Size ii=0; ii < best_state_on_node.size(); ii++ ) {
					if ( best_state_on_node[ii] != 0 )
					TR << rotamer_sets()->rotamer_set_for_moltenresidue( ii+1 )->rotamer( best_state_on_node[ii] )->name1();
					else
					TR << '-';
					}
					TR << ", current_energy: " << currentenergy << ", best_energy: " << bestenergy() << std::endl;
					*/

					bestenergy() = currentenergy;
				}

				if ( record_annealer_trajectory_ ) {
					annealer_trajectory << moltenres_id << " " << rotamer_state_on_moltenres << " A\n";
				}

			} else if ( record_annealer_trajectory_ ) {
				annealer_trajectory << moltenres_id << " " << rotamer_state_on_moltenres << " R\n";
			}
			//else {
			// std::cout << " rejected\n";
			//}


			loopenergy(nn) = currentenergy;
			float const temperature = get_temperature();

			if ( calc_rot_freq() && ( temperature <= calc_freq_temp ) ) {
				++nsteps;
				for ( int ii = 1; ii <= nmoltenres; ++ii ) {
					int iistate = state_on_node(ii);
					if ( iistate != 0 ) {
						++nsteps_for_rot( rotamer_sets()->moltenres_rotid_2_rotid(ii, iistate) );
					}
				}
			}

		} // end of inneriteration loop
	} //end of outeriteration loop

#ifndef NDEBUG
	TR << "pack_rotamers run final: ";
	for ( Size ii=0; ii < best_state_on_node.size(); ii++ ) {
		if ( best_state_on_node[ii] != 0 ) {
			TR << rotamer_sets()->rotamer_set_for_moltenresidue( ii+1 )->rotamer( best_state_on_node[ii] )->name1();
		} else {
			TR << '-';
		}
	}
	TR << ", best_energy: " << bestenergy() << std::endl;
#endif

	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In FixbbSimAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
		std::cerr << "Critical error -- assignment and energy of assignment meaningless" << std::endl;

		FArray1D_int nstates_for_moltenres( rotamer_sets()->nmoltenres(), 0 );
		for ( uint ii = 0; ii < num_rots_to_pack(); ++ii ) {
			++nstates_for_moltenres( rotamer_sets()->res_for_rotamer( rot_to_pack()[ ii ] ) );
		}

		for ( uint ii = 1; ii <= rotamer_sets()->nmoltenres(); ++ii ) {
			if ( best_state_on_node( ii ) == 0 ) {
				std::cout << "Molten res " << ii << " (residue " << rotamer_sets()->moltenres_2_resid( ii );
				std::cout << " ) assigned state 0 despite having " << nstates_for_moltenres( ii ) << " states to choose from" << std::endl;
			}
		}
		debug_assert( ! ig_->any_vertex_state_unassigned() );
		utility_exit();


	}

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( int ii = 1; ii <= nmoltenres; ++ii ) {
		int const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}

	//std::cout << "Annealing ends" << std::endl;
}

void FixbbSimAnnealer::record_annealer_trajectory( bool setting ) {
	record_annealer_trajectory_ = setting;
}

void FixbbSimAnnealer::trajectory_file_name( std::string const & setting ) {
	trajectory_file_name_ = setting;
}


}//end namespace annealer
}//end namespace pack
}//end namespace core
