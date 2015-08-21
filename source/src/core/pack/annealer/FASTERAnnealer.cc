// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FASTERAnnealer.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/annealer/FASTERAnnealer.hh>

// Package headers
#include <core/pack/annealer/FixbbSimAnnealer.hh>
#include <core/pack/interaction_graph/FASTERInteractionGraph.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1A.hh>

// C++ headers
#include <iostream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace annealer {


FASTERAnnealer::FASTERAnnealer(
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::FASTERInteractionGraphOP ig,
	FixbbRotamerSetsCOP rotamer_sets,
	ObjexxFCL::FArray1_int & current_rot_index,
	bool calc_rot_freq,
	ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
):
	parent(
	ig->get_num_total_states(),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	rotamer_sets,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	),
	ig_(ig),
	rotamer_sets_( rotamer_sets ),
	num_nodes_( ig_->get_num_nodes() ),
	recent_network_state_history_( num_nodes_, recent_history_size_, 0 ),
	recent_history_hash_values_( recent_history_size_, 0 ),
	recent_history_hash_count_( hash_size_, 0 ),
	recent_history_head_( 0 ),
	curr_in_recent_history_( 0 ),
	netstate_duplicated_( false ),
	progress_through_sBR_( -1 ),
	sBR_rotamers_( ig_->get_num_total_states() ),
	ciBR_only_( false ),
	num_sa_trajectories_( 5 ),
	sa_inner_iterations_length_scale_( 0.05 ),
	sBR_limit_( -1 )
	//rot_2_moltenres_( 2,ig_->get_num_total_states(), 0)
{
	for ( unsigned int ii =1; ii <= sBR_rotamers_.size(); ++ii ) {
		sBR_rotamers_[ ii-1 ] = ii;
	}
}

FASTERAnnealer::~FASTERAnnealer() {}

void
FASTERAnnealer::run( )
{
	if ( ig_->get_num_total_states() == 0 ) {
		return;
	}

	if ( ciBR_only_ ) {
		ig_->prepare_for_FASTER();
		ig_->blanket_assign_state_0();
		ig_->assign_BMEC();
		if ( ig_->any_vertex_state_unassigned() ) {
			std::cout << "Failed to assign BMEC; some vertices have 0 states" << std::endl;
		}
		//ciBR();
		sBR();
		run_quench_cycles();
		bestenergy() = ig_->get_energy_current_state_assignment();
		finalize_output();
		return;
	}

	core::PackerEnergy best_energy( 0.0 );
	ObjexxFCL::FArray1D_int best_network_state( num_nodes_, 0 );
	ObjexxFCL::FArray1D_int dummy( num_nodes_, 0 );

	//clock_t starttime = clock();

	for ( int ii = 1; ii <= num_sa_trajectories_; ++ii ) {
		{ // scope
			FixbbSimAnnealer fixbb_annealer(
				bestrotamer_at_seqpos(),
				bestenergy(),
				false,
				ig_,
				rotamer_sets_,
				dummy,
				calc_rot_freq(),
				rot_freq()
			);
			fixbb_annealer.scale_inneriterations( sa_inner_iterations_length_scale_ );
			fixbb_annealer.set_assign_state_to_all_nodes_immediately( true );
			fixbb_annealer.run();
		}

		if ( ii == 1 || best_energy > ig_->get_energy_current_state_assignment() ) {
			ig_->get_current_network_state( best_network_state );
			//std::cout << "Best network state:";
			//for ( Size jj = 1; jj <= best_network_state.size(); ++jj ) {
			// std::cout << " " << best_network_state( jj );
			//}
			//std::cout << std::endl;
			best_energy = ig_->get_energy_current_state_assignment();
		}

		ig_->prepare_for_FASTER();
		/*core::PackerEnergy energy = */ig_->get_energy_current_state_assignment();
		//std::cout << "Energy following quick-and-dirty sim annealing: " << energy << std::endl;
		//std::cout << "FASTER::run -- bestenergy: " << best_energy << std::endl;

		sBR();
		core::PackerEnergy energy = ig_->get_energy_current_state_assignment();
		if ( energy < best_energy ) {
			//std::cout << "Found a better energy after completing sBR on rotamer assignment from Fixbb Sim Annealer: " << best_energy << std::endl;
			best_energy = energy;
			ig_->get_current_network_state( best_network_state );
			//for (int jj = 1; jj <= num_nodes_; ++jj)
			//{
			// if ( best_network_state( jj ) < 100 ) std::cout << " ";
			// if ( best_network_state( jj ) < 10  ) std::cout << " ";
			// std::cout << best_network_state( jj ) << " ";
			// if ( jj == 20 ) std::cout << std::endl;
			//}
			//std::cout << std::endl;
		}
	}

	//for (int jj = 1; jj <= num_nodes_; ++jj)
	//{
	// std::cout << "best network state " << std::endl;
	// if ( best_network_state( jj ) < 100 ) std::cout << " ";
	// if ( best_network_state( jj ) < 10  ) std::cout << " ";
	// std::cout << best_network_state( jj ) << " ";
	// if ( jj == 20 ) std::cout << std::endl;
	//}
	//std::cout << std::endl;

	ig_->set_network_state( best_network_state );
	bestenergy() = ig_->get_energy_current_state_assignment();
	//std::cout << "FASTER: bestenergy() " << bestenergy() << std::endl;

	//clock_t stoptime = clock();
	//std::cout << "FASTER completed in " << ((double) stoptime-starttime )/CLOCKS_PER_SEC << " seconds." << std::endl;

	finalize_output();
}


void
FASTERAnnealer::iBR()
{
	reset_recent_network_state_history();
	ObjexxFCL::FArray1D_int current_network_state( num_nodes_, 0 );

	while ( ! stuck_in_network_state_loop() )
			{
		ig_->relax_in_current_context();
		ig_->commit_relaxation();

		//core::PackerEnergy total_energy = ig_->get_energy_current_state_assignment();
		//std::cout << "iBR: total_energy = " << total_energy << std::endl;
		ig_->get_current_network_state( current_network_state );
		note_current_network_state( current_network_state );
	}

}

void
FASTERAnnealer::trySeveral_ciBRs()
{

	ObjexxFCL::FArray1D_int iBRstate( ig_->get_num_nodes() );
	ig_->get_current_network_state( iBRstate );

	ObjexxFCL::FArray1D_int best_ciBRs_state( ig_->get_num_nodes() );
	core::PackerEnergy best_energy = 0;

	for ( int ii = 1; ii < 10; ++ii ) {
		ig_->set_network_state( iBRstate );
		ciBR();
		core::PackerEnergy ciBR_energy = ig_->get_energy_current_state_assignment();
		if ( ii == 1 || ciBR_energy < best_energy ) {
			ig_->get_current_network_state( best_ciBRs_state );
			best_energy = ciBR_energy;
		}
	}
	//std::cout << "ciBR finished with best energy: " << best_energy << std::endl;
	ig_->set_network_state( best_ciBRs_state );
}

void
FASTERAnnealer::ciBR()
{
	Real const DESMET_ciBR_ACCEPT_RATE =  0.8;
	int const limit = 2000; // ??

	reset_recent_network_state_history();
	ObjexxFCL::FArray1D_int current_network_state( num_nodes_, 0 );
	ObjexxFCL::FArray1D_int best_network_state( num_nodes_, 0 );
	core::PackerEnergy ciBR_best_energy( 0 );

	int count = 0;
	while ( ! stuck_in_network_state_loop() ) {
		++count;
		if ( count == limit ) break;

		ig_->relax_in_current_context();
		ig_->probabilistically_commit_relaxation( DESMET_ciBR_ACCEPT_RATE );
		//core::PackerEnergy total_energy = ig_->get_energy_current_state_assignment();
		//std::cout << "ciBR: total_energy = " << total_energy << std::endl;

		ig_->get_current_network_state( current_network_state );

		note_current_network_state( current_network_state );
		if ( count == 1 || ig_->get_energy_current_state_assignment() < ciBR_best_energy ) {
			ciBR_best_energy = ig_->get_energy_current_state_assignment();
			best_network_state = current_network_state;
		}
	}
	ig_->set_network_state( best_network_state );
}

void
FASTERAnnealer::sBR()
{
	shuffle_sBR_rotamers();
	int numAttemptsSinceLastCommit = 0;

	int num_rotamers = static_cast< int >  (sBR_rotamers_.size());

	core::PackerEnergy last_energy = ig_->get_energy_current_state_assignment();
	int count_sBR = 0;
	while ( numAttemptsSinceLastCommit != num_rotamers ) {
		++numAttemptsSinceLastCommit;

		int ran_rotamer = pick_a_rotamer_for_sBR();
		int const node = rotamer_sets_->moltenres_for_rotamer( ran_rotamer );
		int const ran_rotamer_on_node = rotamer_sets_->rotid_on_moltenresidue( ran_rotamer );

		if ( ran_rotamer_on_node == ig_->get_current_state_for_node( node ) ) continue;
		++count_sBR;

		core::PackerEnergy deltaE = ig_->perturb_sBR_and_relax( node, ran_rotamer_on_node );

		//core::PackerEnergy alt_total_energy = ig_->get_energy_following_relaxation();
		//std::cout << "last: " << last_energy << " deltaE: " << deltaE << " last+deltaE: ";
		//std::cout << last_energy + deltaE << " ig's version: " << alt_total_energy << std::endl;


		if ( deltaE < 0.001 ) {
			core::PackerEnergy total_energy = ig_->get_energy_following_relaxation();

			if ( total_energy < last_energy ) {
				//if (deltaE > 0 ) std::cout << "Numerical instability in deltaE calculation?" << deltaE << " " << total_energy << std::endl;

				last_energy = total_energy;
				ig_->commit_relaxation();
				numAttemptsSinceLastCommit = 0;
				shuffle_sBR_rotamers();
				//std::cout << "sBR commit relaxation: total_energy = " << last_energy << std::endl;
			} else {
				ig_->reject_perturbation();
			}

		} else {
			ig_->reject_perturbation();
		}
		if ( sBR_limit_ != -1 && count_sBR >= sBR_limit_ ) {
			break;
		}

	}
	//std::cout << "sBR finished: best_energy = " << last_energy << " in " << count_sBR << " steps; nrot= " << num_rotamers << std::endl;

}

void
FASTERAnnealer::dBR()
{
	core::PackerEnergy last_energy = ig_->get_energy_current_state_assignment();
	for ( int ii = 1; ii < 10000; ++ii ) {
		int node1 = ((int) ( numeric::random::rg().uniform() * num_nodes_)) + 1;
		int node2 = ig_->get_random_neighbor_for_node( node1 );

		if ( node2 == 0 ) continue;

		int rot1 = pick_rotamer_for_node( node1 );
		int rot2 = pick_rotamer_for_node( node2 );

		core::PackerEnergy deltaE  = ig_->perturb_dBR_and_relax( node1, rot1, node2, rot2 );

		if ( deltaE < 0.001 ) {
			core::PackerEnergy total_energy = ig_->get_energy_following_relaxation();
			if ( total_energy < last_energy ) {
				last_energy = total_energy;
				ig_->commit_relaxation();
				//std::cout << "dBR commit relaxation: total_energy = " << last_energy << std::endl;
			} else {
				ig_->reject_perturbation();
			}
		} else {
			ig_->reject_perturbation();
		}
	}
	//std::cout << "dBR finished: best_energy = " << last_energy << std::endl;
}

void FASTERAnnealer::set_ciBR_only( bool setting )
{
	ciBR_only_ = setting;
}

void FASTERAnnealer::set_num_sa_trajectories( int setting )
{
	num_sa_trajectories_ = setting;
}

void FASTERAnnealer::set_sa_length_scale( Real setting )
{
	sa_inner_iterations_length_scale_ = setting;
}


void
FASTERAnnealer::finalize_output()
{

	ObjexxFCL::FArray1D_int best_state_on_node( num_nodes_, 0 );
	ig_->get_current_network_state( best_state_on_node );

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( Size ii = 1; ii <= rotamer_sets_->nmoltenres(); ++ii ) {
		int const iiresid = rotamer_sets_->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets_->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}

}

void
FASTERAnnealer::reset_recent_network_state_history()
{
	recent_network_state_history_ = 0;
	curr_in_recent_history_ = 0;
	recent_history_head_ = 0;
	recent_history_hash_values_ = 0;
	recent_history_hash_count_ = 0;
	netstate_duplicated_ = false;
}

void
FASTERAnnealer::note_current_network_state( ObjexxFCL::FArray1_int const & netstate )
{
	if ( curr_in_recent_history_ < recent_history_size_ ) {
		++recent_history_head_;
		++curr_in_recent_history_;

	} else {
		++recent_history_head_;
		if ( recent_history_head_ > recent_history_size_ ) {
			recent_history_head_ = 1;
		}
		--recent_history_hash_count_( recent_history_hash_values_( recent_history_head_ ) );

	}

	ObjexxFCL::FArray1A_int history_line( recent_network_state_history_( 1, recent_history_head_ ), num_nodes_ );
	history_line = netstate;

	int hash = hash_recent_history( recent_history_head_ );
	int num_same_hash = recent_history_hash_count_( hash );
	if ( num_same_hash != 0 ) {
		//possible duplicate: stupid, slow method -- look at them all
		for ( int ii = 1; ii <= curr_in_recent_history_; ++ii ) {
			if ( recent_history_head_ == ii ) continue;
			if ( hash == recent_history_hash_values_( ii ) ) {
				--num_same_hash;
				bool same = true;
				for ( int jj = 1; jj <= num_nodes_; ++jj ) {
					if ( recent_network_state_history_( jj, recent_history_head_ ) !=
							recent_network_state_history_( jj, ii ) ) {
						same = false;
						break;
					}
				}
				if ( same ) {
					netstate_duplicated_ = true;
					break;
				}

				if ( num_same_hash == 0 ) {
					break;
				}
			}
		}
	}

	recent_history_hash_values_( recent_history_head_ ) = hash;
	++recent_history_hash_count_( hash );

}

bool
FASTERAnnealer::stuck_in_network_state_loop()
{
	return netstate_duplicated_;
}

void
FASTERAnnealer::run_quench_cycles()
{
	ObjexxFCL::FArray1D_int state_on_node( ig_->get_num_nodes(), 0 );
	ig_->get_current_network_state( state_on_node );
	int count_since_last_commit( 0 );
	int count_total( 0 );
	int num_rotamers = static_cast< int >  (sBR_rotamers_.size());

	core::PackerEnergy deltaE, dummy, best_energy( ig_->get_energy_current_state_assignment() );

	set_to_quench();


	while ( count_since_last_commit < num_rotamers ) {
		++count_since_last_commit;
		++count_total;

		int ran_rotamer = pick_a_rotamer( count_total );
		int const node = rotamer_sets_->moltenres_for_rotamer( ran_rotamer );
		int const ran_rotamer_on_node = rotamer_sets_->rotid_on_moltenresidue( ran_rotamer );

		if ( state_on_node( node ) == ran_rotamer_on_node ) continue;

		ig_->consider_substitution( node, ran_rotamer_on_node, deltaE, dummy );
		if ( state_on_node( node ) == 0 || deltaE < 0 ) {
			core::PackerEnergy totalE = ig_->commit_considered_substitution();
			state_on_node( node ) = ran_rotamer_on_node;
			if ( totalE + 1e-4 < best_energy ) {
				count_since_last_commit = 0;
				best_energy = totalE;
			}
		}

		if ( count_since_last_commit > 2*num_rotamers ) {
			break;
		}
	}

	set_not_to_quench();

}


int
FASTERAnnealer::pick_rotamer_for_node( int node )
{
	int num_states_for_node = ig_->get_num_states_for_node( node );
	int rand_state = ((int) ( numeric::random::rg().uniform() * num_states_for_node ) ) + 1;
	return rand_state;
}

int
FASTERAnnealer::pick_a_rotamer_for_sBR()
{
	++progress_through_sBR_;
	return sBR_rotamers_[ progress_through_sBR_ ];
}

void
FASTERAnnealer::shuffle_sBR_rotamers()
{
	progress_through_sBR_ = -1;
	numeric::random::random_permutation( sBR_rotamers_.begin(), sBR_rotamers_.end(), numeric::random::rg() );
}

int
FASTERAnnealer::hash_recent_history( int history_index )
{
	int hash = 0;
	for ( int ii = 1; ii <= num_nodes_; ++ii ) {
		hash += ii * recent_network_state_history_( ii, history_index );
		hash = hash % hash_size_;
	}
	return hash;
}

}
}
}


