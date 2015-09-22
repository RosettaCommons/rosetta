// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/MultiCoolAnnealer.cc
/// @brief  Multiple low-temperature cooling cycles annealer class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <core/pack/annealer/MultiCoolAnnealer.hh>

/// Package headers
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/task/PackerTask.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/FArray1A.hh>

/// Utility headers
#include <utility/exit.hh>

// Numeric headers
#include <numeric/random/random.hh>

/// C++ headers
#include <iostream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#ifdef WIN32
#include <ctime>
#endif


namespace core {
namespace pack {
namespace annealer {

static THREAD_LOCAL basic::Tracer TR( "core.pack.annealer.MultiCoolAnnealer" );

using namespace ObjexxFCL;
using namespace pack::interaction_graph;
using namespace pack::rotamer_set;
using namespace pack::task;

core::PackerEnergy const MultiCoolAnnealer::uninitialized_energy( 1234 );

/// @brief constructor
MultiCoolAnnealer::MultiCoolAnnealer(
	PackerTaskCOP task,
	utility::vector0< int > & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	InteractionGraphBaseOP ig,
	FixbbRotamerSetsCOP p_rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq
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
	),
	ig_(ig),
	nsteps_for_rot_( p_rotamer_set->nrotamers() , 0 ),
	nsteps_( 0 ),
	top_to_keep( task->multi_cool_annealer_history_size() ),
	top_netstates_( ig_->get_num_nodes(), top_to_keep, 0 ),
	energy_top_( top_to_keep, uninitialized_energy ),
	worst_top_energy_( uninitialized_energy ),
	which_netstate_worst_top_( 1 ),
	num_top_kept_( 0 )
{
}

MultiCoolAnnealer::MultiCoolAnnealer(
	PackerTaskCOP task,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	InteractionGraphBaseOP ig,
	FixbbRotamerSetsCOP p_rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq
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
	),
	ig_(ig),
	nsteps_for_rot_( num_rots_to_pack(), 0 ),
	nsteps_( 0 ),
	top_to_keep( task->multi_cool_annealer_history_size() ),
	top_netstates_( ig_->get_num_nodes(), top_to_keep, 0 ),
	energy_top_( top_to_keep, uninitialized_energy ),
	worst_top_energy_( uninitialized_energy ),
	which_netstate_worst_top_( 1 ),
	num_top_kept_( 0 )
{
}

/// @brief virtual destructor
MultiCoolAnnealer::~MultiCoolAnnealer()
{}


/// @brief sim_annealing for fixed backbone design mode
void MultiCoolAnnealer::run()
{

	//--------------------------------------------------------------------
	//internal variables

	//FArray1D_int list( rotamer_sets()->nrotamers() );
	FArray1D_int state_on_node( rotamer_sets()->nmoltenres(),0 );
	FArray1D_int best_state_on_node( rotamer_sets()->nmoltenres(),0 );


	//bk variables for calculating rotamer frequencies during simulation

	FArray1D_int network_state_to_restore( rotamer_sets()->nmoltenres() );
	FArray1D< core::PackerEnergy > temp_top_generated_at( top_to_keep, 0.0f );
	FArray1D< core::PackerEnergy > second_round_cooling_finalE( top_to_keep, 0.0f );
	FArray2D_int hamming_distance( top_to_keep, top_to_keep, 0 );
	FArray1D_int hamming_start_to_stop( top_to_keep, 0 );

	//unsigned int count_considered_rotamer_subs = 0;

	//--------------------------------------------------------------------
	//initialize variables


	ig_->prepare_for_simulated_annealing();

	TR << "IG after prepare_for_simulated_annealing: " << ig_->getTotalMemoryUsage() << " bytes" << std::endl;

	ig_->blanket_assign_state_0();

	//--------------------------------------------------------------------
	if ( num_rots_to_pack() == 0 ) return;


	/// Start off in some random state assignment to avoid unassigned-state problems that
	/// are sometimes seen in optH
	for ( int ii = 1; ii <= ig_->get_num_nodes(); ++ii ) {
		float dummy_delta( 0.0 ), dummy_prev_energy( 0.0 );
		if ( ig_->get_num_states_for_node( ii ) == 1 ) {
			ig_->consider_substitution( ii, 1, dummy_delta, dummy_prev_energy );
			ig_->commit_considered_substitution();
			state_on_node( ii ) = 1;
		} else {
			int randstate = rotamer_sets()->rotid_on_moltenresidue( pick_a_rotamer_for_node( ii ));
			ig_->consider_substitution( ii, randstate, dummy_delta, dummy_prev_energy );
			ig_->commit_considered_substitution();
			state_on_node( ii ) = randstate;
		}
	}
	bestenergy() = ig_->get_energy_current_state_assignment();
	best_state_on_node = state_on_node;


	setup_iterations();

	int outeriterations = get_outeriterations() * 6;
	//core::PackerEnergy last_temperature = 10;
	//outer loop
	set_temperature( 10 );
	set_lowtemp( 0.2 );
	for ( int nn = 1; nn <= outeriterations; ++nn ) {

		if ( nn % 6 == 1 && nn != 1 ) {
			cool();
		}

		int inneriterations = get_inneriterations();
		inneriterations /= 5;

		//std::cerr << "inneriterations: " << inneriterations << std::endl;

		if ( nn % 2 == 0 ) {
			//last_temperature = get_temperature();
			network_state_to_restore = state_on_node;

			bestenergy() = ig_->get_energy_current_state_assignment();
			best_state_on_node = state_on_node;
			run_quench( state_on_node, best_state_on_node, bestenergy(), inneriterations );
			store_top_energy( best_state_on_node, bestenergy() );

			state_on_node = network_state_to_restore;
			ig_->set_network_state( state_on_node );
		} else {
			run_constant_temp_rotamer_substitutions(
				state_on_node,
				best_state_on_node,
				bestenergy(),
				inneriterations);
			store_top_energy( best_state_on_node, bestenergy() );
			//std::cerr << "Temperature: " << get_temperature() << " bestenergy so far: " << bestenergy() << std::endl;
		}

		//std::cerr << "Considered " << nsteps_ << " rotamer substitutions so far." << std::endl;
	} //end of outeriteration loop

	{  // scope

		core::PackerEnergy best_of_best = energy_top_( 1 );
		int which_best_of_best = 1;
		for ( Size ii = 2; ii <= top_to_keep; ++ii ) {
			if ( best_of_best > energy_top_( ii ) && energy_top_( ii ) != uninitialized_energy ) {
				best_of_best = energy_top_( ii );
				which_best_of_best = ii;
			}
		}
		FArray1A_int best_of_best_state( top_netstates_( 1, which_best_of_best), rotamer_sets()->nmoltenres() );
		best_state_on_node = best_of_best_state;
		bestenergy() = best_of_best;
	}

	set_lowtemp( 0.05 );
	for ( Size ii = 1; ii <= top_to_keep; ++ii ) {
		//clock_t starttime = clock();
		//std::cout << "MultiCoolAnnealer: starting low temp annealing # "<< ii  << " of " << top_to_keep << std::endl;
		//std::cout << " with network state: " << std::endl;
		//for ( int jj = 1; jj <= rotamer_sets()->nmoltenres(); ++jj )
		//{
		// int jjstate = top_netstates_(jj, ii );
		// if ( jjstate < 1000 ) std::cerr << " ";
		// if ( jjstate < 100  ) std::cerr << " ";
		// if ( jjstate < 10   ) std::cerr << " ";
		// std::cout << jjstate << " ";
		// if ( jj % 10 == 0 ) std::cerr << std::endl;
		//}
		//std::cout << std::endl << "With energy: " << energy_top_( ii ) << std::endl;

		if ( energy_top_( ii ) == uninitialized_energy ) continue;

		FArray1A_int start_state( top_netstates_( 1, ii ), rotamer_sets()->nmoltenres() );
		state_on_node = start_state;
		ig_->set_network_state( state_on_node );
		core::PackerEnergy best_energy_this_starting_point = energy_top_( ii );
		FArray1D_int best_state_on_node_this_starting_point( rotamer_sets()->nmoltenres() );
		best_state_on_node_this_starting_point = state_on_node;

		set_temperature( 0.25 );
		outeriterations = 6;
		int inneriterations = get_inneriterations();

		for ( int nn = 1; nn <= outeriterations; ++nn ) {
			if ( nn != 1 ) cool();

			//std::cerr << "MultiCoolAnnealer: temperature = " << get_temperature() << " currenergy: ";
			//std::cerr << ig_->get_energy_current_state_assignment() << std::endl;

			run_constant_temp_rotamer_substitutions(
				state_on_node,
				best_state_on_node_this_starting_point,
				best_energy_this_starting_point,
				inneriterations);

			FArray1D_int state_to_restore( state_on_node );
			//core::PackerEnergy bestE_before_quench( bestenergy() );
			run_quench(
				state_on_node,
				best_state_on_node,
				bestenergy(),
				inneriterations / 2);
			//if ( bestE_before_quench > bestenergy() )
			//{
			// std::cerr << "Quench inside cooling run # " << ii << " found lower energy: ";
			// std::cerr << bestenergy() << std::endl;
			//}
			state_on_node = state_to_restore;
			ig_->set_network_state( state_on_node );

		} // end outer iterations

		if ( std::abs(best_energy_this_starting_point - energy_top_( ii ) ) < 0.001 ) {
			bool same_as_start = true;
			for ( Size jj = 1; jj <= rotamer_sets()->nmoltenres(); ++jj ) {
				if ( best_state_on_node_this_starting_point( jj ) != top_netstates_( jj, ii ) ) {
					same_as_start = false;
					break;
				}
			}

			if ( same_as_start ) {
				//clock_t stoptime = clock();
				//std::cout << "Round " << ii << " cooling completed in " << ((double) stoptime-starttime)/CLOCKS_PER_SEC << std::endl;
				continue;
			}
		}

		state_on_node = best_state_on_node_this_starting_point;
		ig_->set_network_state( state_on_node );
		run_quench( state_on_node,
			best_state_on_node_this_starting_point,
			best_energy_this_starting_point,
			inneriterations );

		second_round_cooling_finalE( ii ) = best_energy_this_starting_point;
		for ( Size jj = 1; jj <= rotamer_sets()->nmoltenres(); ++jj ) {
			if ( best_state_on_node_this_starting_point( jj ) != top_netstates_(jj, ii ) ) {
				++hamming_start_to_stop(ii);
			}
		}

		if ( bestenergy() > best_energy_this_starting_point ) {
			//std::cerr << "BEST FOUND" << std::endl;
			bestenergy() = best_energy_this_starting_point;
			best_state_on_node = best_state_on_node_this_starting_point;
		}

		//std::cerr << "MultiCoolAnnealer: final low temp annealing # "<< ii << " with network state: " << std::endl;
		//for ( int jj = 1; jj <= rotamer_sets()->nmoltenres(); ++jj )
		//{
		// int jjstate = state_on_node(jj);
		// if ( jjstate < 1000 ) std::cerr << " ";
		// if ( jjstate < 100  ) std::cerr << " ";
		// if ( jjstate < 10   ) std::cerr << " ";
		// std::cerr << jjstate << " ";
		// if ( jj % 10 == 0 ) std::cerr << std::endl;
		//}
		//std::cerr << std::endl << "With energy: " << ig_->get_energy_current_state_assignment();
		//std::cerr << " and bestenergy()" << bestenergy() << std::endl;

		//clock_t stoptime = clock();
		//std::cout << "Round " << ii << " cooling completed in " << ((double) stoptime-starttime)/CLOCKS_PER_SEC << std::endl;
	} // end top_states

	ig_->set_network_state( best_state_on_node );

	//std::cerr << "MultiCoolAnnealer finished with an energy of " << bestenergy() << " after considering " << nsteps_ << " rotamer substitutions. " << std::endl;

	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In MultiCoolAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
		std::cerr << "Critical error -- assignment and energy of assignment meaningless" << std::endl;

		FArray1D_int nstates_for_moltenres( rotamer_sets()->nmoltenres(), 0 );
		for ( Size ii = 0; ii < num_rots_to_pack(); ++ii ) {
			++nstates_for_moltenres( rotamer_sets()->res_for_rotamer( rot_to_pack()[ ii ] ) );
		}

		for ( Size ii = 1, iie = rotamer_sets()->nmoltenres(); ii <= iie; ++ii ) {
			if ( best_state_on_node( ii ) == 0 ) {
				std::cerr << "Molten res " << ii << " (residue " << rotamer_sets()->moltenres_2_resid( ii );
				std::cerr << " ) assigned state 0 despite having " << nstates_for_moltenres( ii ) << " states to choose from" << std::endl;
			}
		}
		std::cout << "num_top_kept_: " << num_top_kept_ << std::endl;
		debug_assert( ! ig_->any_vertex_state_unassigned() );
		utility_exit();
	}

	//for (int ii = 1; ii <= top_to_keep; ++ii)
	//{
	// for (int jj = 1; jj <= ii; ++jj )
	// {
	//  std::cerr << "    ";
	// }
	//
	// for (int jj = ii + 1; jj <= top_to_keep; ++jj)
	// {
	//  int hamming = 0;
	//  for (int kk = 1; kk <= rotamer_sets()->nmoltenres(); ++kk)
	//  {
	//   if ( top_netstates_( kk, ii ) != top_netstates_( kk, jj ) )
	//   {
	//    ++hamming;
	//   }
	//  }
	//  hamming_distance( jj, ii ) = hamming;
	//  if ( hamming < 100  ) std::cerr << " ";
	//  if ( hamming < 10   ) std::cerr << " ";
	//  std::cerr << hamming << " ";
	// }
	//
	// std::cerr << energy_top_(ii) << " " << hamming_start_to_stop(ii) << " "<< second_round_cooling_finalE( ii ) << std::endl;
	//
	//}

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( Size ii = 1; ii <= rotamer_sets()->nmoltenres(); ++ii ) {
		int iiresid = rotamer_sets()->moltenres_2_resid(ii);
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii) );
	}
}

void
MultiCoolAnnealer::cool()
{
	core::PackerEnergy temperature = get_temperature();
	temperature = ( temperature - get_lowtemp() ) * std::exp( -1. ) + get_lowtemp();
	set_temperature( temperature );
}

void
MultiCoolAnnealer::run_quench
(
	FArray1D_int & state_on_node,
	FArray1D_int & best_state_on_node,
	core::PackerEnergy & best_energy,
	int num_cycles
)
{
	set_to_quench();
	run_constant_temp_rotamer_substitutions(
		state_on_node, best_state_on_node,
		best_energy, num_cycles );
	set_not_to_quench();
}

void
MultiCoolAnnealer::run_constant_temp_rotamer_substitutions(
	FArray1D_int & state_on_node,
	FArray1D_int & best_state_on_node,
	core::PackerEnergy & best_energy,
	int num_cycles
)
{

	core::PackerEnergy currentenergy = ig_->get_energy_current_state_assignment();
	int substitutions_without_a_commit = 0;

	core::PackerEnergy threshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
	if ( quench() ) { threshold_for_deltaE_inaccuracy = 0; }
	ig_->set_errorfull_deltaE_threshold( threshold_for_deltaE_inaccuracy );

	for ( int n = 1; n <= num_cycles; ++n ) {

		int ranrotamer = pick_a_rotamer( n );
		if ( ranrotamer == -1 ) continue;

		int rotamer_seqpos = rotamer_sets()->res_for_rotamer(ranrotamer);
		int moltenres_id = rotamer_sets()->resid_2_moltenres(rotamer_seqpos);
		int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue(ranrotamer);
		int prevrotamer_state = state_on_node(moltenres_id);

		if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

		core::PackerEnergy delta_energy, previous_energy_for_node;
		++nsteps_;
		ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
			delta_energy, previous_energy_for_node);

		//bk keep new rotamer if it is lower in energy or accept it at some
		//bk probability if it is higher in energy, if it is the first
		//bk rotamer to be tried at this position automatically accept it.
		if ( (prevrotamer_state == 0)||pass_metropolis(previous_energy_for_node,delta_energy) ) {
			substitutions_without_a_commit = 0;
			currentenergy = ig_->commit_considered_substitution();
			state_on_node(moltenres_id) = rotamer_state_on_moltenres;
			if ( (prevrotamer_state == 0)||(currentenergy < best_energy) ) {
				if ( ! ig_->any_vertex_state_unassigned() ) {
					best_energy = currentenergy;
					best_state_on_node = state_on_node;
				}
			}
		} else {
			++substitutions_without_a_commit;
			if ( quench() && (Size) substitutions_without_a_commit == 2 * num_rots_to_pack() ) {
				//visited all rotamers without finding any that lowered the energy
				return;
			}
		}// end Metropolis criteria

		if ( calc_rot_freq() && ( get_temperature() <= calc_freq_temp ) ) {
			for ( Size ii = 1; ii <= rotamer_sets()->nmoltenres(); ++ii ) {
				int iistate = state_on_node(ii);
				if ( iistate != 0 ) {
					++nsteps_for_rot_( rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii) ) );
				}
			}
		}

	} //end inner iterations
}


void MultiCoolAnnealer::store_top_energy(
	FArray1D_int const & state_on_node,
	core::PackerEnergy energy
)
{
	if ( worst_top_energy_ != uninitialized_energy && energy > worst_top_energy_ ) return;
	for ( Size ii = 1; ii <= rotamer_sets()->nmoltenres(); ++ii ) {
		if ( state_on_node(ii) == 0 ) return;
	}

	for ( Size ii = 1; ii <= top_to_keep; ++ii ) {
		if ( std::abs( energy - energy_top_( ii )) < 0.0001 ) {
			bool repeat = true;
			for ( Size jj = 1; jj <= rotamer_sets()->nmoltenres(); ++jj ) {
				if ( state_on_node( jj ) != top_netstates_( jj, ii ) ) {
					repeat = false;
					break;
				}
			}
			if ( repeat ) {
				//std::cerr << "MultiCoolAnnealer: already found this network state" << std::endl;
				return;
			}
		}
	}

	//std::cerr << "MultiCoolAnnealer:  replacing worst_best: " <<
	// which_netstate_worst_top_ << " " << worst_top_energy_ << " with " <<
	// energy << std::endl;

	FArray1A_int netstate_to_replace(
		top_netstates_( 1, which_netstate_worst_top_ ), rotamer_sets()->nmoltenres() );

	//keep this netstate
	netstate_to_replace = state_on_node;
	energy_top_( which_netstate_worst_top_ ) = energy;
	++num_top_kept_;

	worst_top_energy_ = energy_top_(1);
	which_netstate_worst_top_ = 1;

	for ( Size ii = 2; ii <= top_to_keep; ++ii ) {
		if ( worst_top_energy_ < energy_top_( ii ) || energy_top_( ii ) == uninitialized_energy ) {
			worst_top_energy_ = energy_top_( ii );
			which_netstate_worst_top_ = ii;
		}
		if ( worst_top_energy_ == uninitialized_energy ) break;
	}
	//std::cerr << "next worst_top : " << worst_top_energy_ << " from #" << which_netstate_worst_top_ << std::endl;

}

} /// namespace annealer
} // namespace pack
} // namespace core
