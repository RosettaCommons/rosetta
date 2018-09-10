// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FixbbPwatSimAnnealer.cc
/// @brief  fixed temperature annealer for packing with statistical PWAT water model
/// @author  Ryan Pavlovicz (rpavlov@uw.edu)

// Unit Headers
#include <core/pack/annealer/FixbbPwatSimAnnealer.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack_basic/RotamerSetsBase.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>

#include <iostream>
#include <fstream>

#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace annealer {

////////////////////////////////////////////////////////////////////////////////
FixbbPwatSimAnnealer::FixbbPwatSimAnnealer(
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
	ig_(ig)
{
	min_dwell_ = 0.02;
}

/// @brief virtual destructor
FixbbPwatSimAnnealer::~FixbbPwatSimAnnealer()
{}

/// @brief sim_annealing for fixed backbone design mode
void FixbbPwatSimAnnealer::run() {
	int const nmoltenres = ig_->get_num_nodes();
	int totalrot = 0;
	for ( int resi = 1; resi <= nmoltenres; resi++ ) {
		totalrot += rotamer_sets()->nrotamers_for_moltenres(resi);
	}

	// the temperature in steady state versus spikes
	core::Real temp_low=1.0, temp_high=100.0;

	FArray1D_int state_on_node( nmoltenres,0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D_float loopenergy(maxouteriterations,0.0);

	//--------------------------------------------------------------------
	//initialize variables
	core::PackerEnergy currentenergy = 0.0;

	ig_->prepare_for_simulated_annealing();
	ig_->blanket_assign_state_0();
	//--------------------------------------------------------------------
	if ( num_rots_to_pack() == 0 ) return;

	setup_iterations();

	int outeriterations = 50; // fd override base class (the meaning of outer iterations is somewhat different...)

	// recording info
	core::Size record_count = 0;
	utility::vector1< utility::vector1<core::Size> > rot_states(nmoltenres);
	for ( int roti=1; roti<=nmoltenres; roti++ ) {
		rot_states[roti].resize(rotamer_sets()->nrotamers_for_moltenres(roti),0);
	}

	// outer loop
	for ( int nn = 1; nn <= outeriterations; ++nn ) {
		if ( quench() ) {
			currentenergy = bestenergy();
			state_on_node = best_state_on_node;
			ig_->set_network_state( state_on_node );
		}

		float treshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( treshold_for_deltaE_inaccuracy );

		// do high temperature spike!
		set_temperature( temp_high );
		int spike_steps = totalrot;

		for ( int n = 1; n <= spike_steps; ++n ) {
			int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue;

			core::conformation::Residue const & curres( *rotamer_sets()->rotamer_for_moltenres(moltenres_id, rotamer_state_on_moltenres) );
			if ( curres.aa() == core::chemical::aa_h2o ) {
				double rand_num = numeric::random::rg().uniform();
				double nrot = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
				if ( rand_num < (nrot/2.0-1)/nrot ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
			}
			core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );
			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
				delta_energy, previous_energy_for_node);
			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					best_state_on_node = state_on_node;
					bestenergy() = currentenergy;
				}
			}
		}

		// back to low temperature
		set_temperature( temp_low );
		int inneriterations = 6*totalrot;

		int accepts = 0;
		for ( int n = 1; n <= inneriterations; ++n ) {  // normal inneriterations
			int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// for waters, set to virtual 50% of the time
			core::conformation::Residue const &curres( *rotamer_sets()->rotamer_for_moltenres(moltenres_id, rotamer_state_on_moltenres) );
			if ( curres.aa() == core::chemical::aa_h2o ) {
				double rand_num = numeric::random::uniform();
				// last rotamer is the virtual state
				double nrot = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
				if ( rand_num < 0.5*(1.0 - 1.0/nrot) ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
			}

			// initializing to zero but should be updated below.
			core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );

			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
				delta_energy, previous_energy_for_node);

			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
				accepts++;
				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					best_state_on_node = state_on_node;
					bestenergy() = currentenergy;
				}
			}

			// accumulate dwell times
			if ( n > totalrot ) {  // allow for totalrot steps of equilibration before recording
				record_count++;
				for ( int roti=1; roti <= nmoltenres; ++roti ) {
					if ( state_on_node(roti) == 0 ) {
						continue;  // this rot hasn't been visited yet (this should not happen much...)
					} else {
						rot_states[roti][state_on_node(roti)]++;
					}
				}
			}

			loopenergy(nn) = currentenergy;
		} // end of inneriteration loop
	} //end of outeriteration loop

	all_rot_.clear();
	for ( Size x = 1; x <= rot_states.size(); ++x ) {
		std::string resname = rotamer_sets()->rotamer_set_for_moltenresidue(x)->rotamer(1)->name(); // better way?
		if ( resname != "PWAT" ) continue;

		for ( Size y = 1; y < rot_states[x].size(); ++y ) { // last res is VRT, ignore
			PointDwell current_rotamer;
			current_rotamer.dwell = (core::Real)rot_states[x][y] / (core::Real)record_count;
			core::conformation::ResidueCOP pwat_rot = rotamer_sets()->rotamer_set_for_moltenresidue(x)->rotamer(y);
			current_rotamer.xyz = pwat_rot->xyz("O");
			if ( current_rotamer.dwell >= min_dwell_ ) {
				all_rot_.push_back(current_rotamer);
			}
		}
	}

	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In FixbbPwatSimAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
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
		assert( ! ig_->any_vertex_state_unassigned() );
		utility_exit();
	}

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( int ii = 1; ii <= nmoltenres; ++ii ) {
		int const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}

	//std::cout << "Annealing ends" << std::endl;
}

}//end namespace annealer
}//end namespace pack
}//end namespace core
