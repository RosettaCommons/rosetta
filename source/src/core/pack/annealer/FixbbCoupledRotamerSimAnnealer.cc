// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FixbbCoupledRotamerSimAnnealer.cc
/// @brief  Packer's standard simulated annealing class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/FixbbCoupledRotamerSimAnnealer.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <numeric/random/random.hh>


#include <iostream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace annealer {

static THREAD_LOCAL basic::Tracer TR( "core.pack.annealer.FixbbCoupledRotamerSimAnnealer", basic::t_info );

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
FixbbCoupledRotamerSimAnnealer::FixbbCoupledRotamerSimAnnealer(
	utility::vector0<int> & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_sets,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq,
	RotamerCouplingsCOP rotamer_couplings
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
	), ig_(ig)
{
	setup_rotamer_couplings( rotamer_couplings );
}

FixbbCoupledRotamerSimAnnealer::FixbbCoupledRotamerSimAnnealer(
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq,
	RotamerCouplingsCOP rotamer_couplings
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
	), ig_(ig)
{
	setup_rotamer_couplings( rotamer_couplings );
}


void
FixbbCoupledRotamerSimAnnealer::setup_rotamer_couplings(
	RotamerCouplingsCOP rotamer_couplings
)
{
	// setup the rotamer couplings
	// for each moltenres
	rotamer_couplings_ = RotamerCouplingsOP( new rotamer_set::RotamerCouplings() );
	rotamer_couplings_->resize( rotamer_sets()->nmoltenres() );

	for ( Size moltenres_id=1; moltenres_id<= rotamer_sets()->nmoltenres(); ++moltenres_id ) {
		uint const       resid( rotamer_sets()->moltenres_2_resid( moltenres_id ) );
		uint const other_resid( (*rotamer_couplings)[ resid ].first );

		if ( other_resid && rotamer_sets()->resid_2_moltenres( other_resid ) ) {
			// both are molten
			(*rotamer_couplings_)[ moltenres_id ].first  = rotamer_sets()->resid_2_moltenres( other_resid );
			(*rotamer_couplings_)[ moltenres_id ].second = (*rotamer_couplings)[ resid ].second;
		} else {
			(*rotamer_couplings_)[ moltenres_id ].first = 0; // signal no pairwise couplings
		}
	}

}


/// @brief virtual destructor
FixbbCoupledRotamerSimAnnealer::~FixbbCoupledRotamerSimAnnealer()
{}


void FixbbCoupledRotamerSimAnnealer::run()
{
	int const nmoltenres = ig_->get_num_nodes();

	FArray1D_int state_on_node( nmoltenres,0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D< core::PackerEnergy > loopenergy(maxouteriterations,0.0);

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

	FArray1D< core::PackerEnergy > previous_nsteps_for_rot( rotamer_sets()->nrotamers(), 0.0);

	int outeriterations = get_outeriterations();

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

		core::PackerEnergy treshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( treshold_for_deltaE_inaccuracy );

		//inner loop
		for ( int n = 1; n <= inneriterations; ++n ) {
			int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int const rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			core::PackerEnergy previous_energy_for_node, delta_energy;

			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
				delta_energy, previous_energy_for_node);

			// specialize to the case of coupled pairs in this first pass implementation

			// assume all couplings are between moltenres -- couplings between fixed and moltenres could
			// have been used to trim the residueset

			int const other_moltenres_id( (*rotamer_couplings_)[ moltenres_id ].first );

			if ( other_moltenres_id >= 1 ) {
				TR.Trace << "moltenres_id " << moltenres_id << " coupled to moltenres_id " << other_moltenres_id << std::endl;
				conformation::ResidueMatcherCOP matcher( (*rotamer_couplings_)[ moltenres_id ].second );

				RotamerSetCOP rotamer_set( rotamer_sets()->rotamer_set_for_moltenresidue( moltenres_id ) );
				conformation::ResidueCOP new_rotamer( rotamer_set->rotamer( rotamer_state_on_moltenres ) );

				// check if new rotamer is compatible with current rotamer
				int const other_prevrotamer_state( state_on_node( other_moltenres_id ) );
				RotamerSetCOP other_rotamer_set( rotamer_sets()->rotamer_set_for_moltenresidue( other_moltenres_id ) );
				conformation::ResidueCOP other_rotamer( other_prevrotamer_state == 0 ? conformation::ResidueCOP(0) :
					other_rotamer_set->rotamer( other_prevrotamer_state ) );

				if ( !( other_rotamer && (*matcher)( *new_rotamer, *other_rotamer ) ) ) {
					// need to try a double substitution!!
					TR.Trace << "Picked rotamer incompatible, trying double substitution" << std::endl;
					int other_rotamer_state(0);
					int const other_nrotamers( other_rotamer_set->num_rotamers() );
					int tries(0);
					while ( true ) {
						++tries;
						if ( tries > 1000 ) TR.Trace << "tries: " << tries << '\n';
						// pick a new rotamer at the other position
						other_rotamer_state = static_cast<int>( other_nrotamers * numeric::random::rg().uniform() + 1 );
						other_rotamer = other_rotamer_set->rotamer( other_rotamer_state );

						if ( (*matcher)(*new_rotamer, *other_rotamer ) ) {
							TR.Trace << "matching rotamer found" << std::endl;
							break;
						}
					}

					// now make the double substitution
					core::PackerEnergy tmp_currentenergy = ig_->commit_considered_substitution();
					{ // debugging
						Real const dev( std::abs( tmp_currentenergy - currentenergy - delta_energy ) );
						if ( dev > 0.01 ) {
							TR << "equal1? " << dev << ' ' << tmp_currentenergy - currentenergy << ' ' << delta_energy <<
								'\n';
						}
					} // scope

					core::PackerEnergy delta_energy2, previous_energy_for_node2;
					ig_->consider_substitution( other_moltenres_id, other_rotamer_state,
						delta_energy2, previous_energy_for_node2 );

					core::PackerEnergy previous_energy_average = ( previous_energy_for_node + previous_energy_for_node2 ) / 2.0;
					core::PackerEnergy delta_energy_average = ( delta_energy + delta_energy2 ) / 2.0;

					if ( prevrotamer_state == 0 || other_prevrotamer_state == 0 ||
							pass_metropolis( previous_energy_average, delta_energy_average ) ) {
						// accept !!!!!!!
						TR.Trace << "accepting coupled rotamer substitution" << std::endl;
						currentenergy = ig_->commit_considered_substitution();

						{ // debugging
							Real const dev( std::abs( currentenergy - tmp_currentenergy - delta_energy2 ) );
							if ( dev > 0.01 ) {
								TR << "equal2? " << dev << ' ' << currentenergy - tmp_currentenergy << ' ' << delta_energy2 <<
									'\n';
							}
						} // scope

						state_on_node(       moltenres_id ) = rotamer_state_on_moltenres;
						state_on_node( other_moltenres_id ) = other_rotamer_state;


						if ( ( prevrotamer_state == 0 ) || ( other_prevrotamer_state == 0 ) || ( currentenergy < bestenergy() ) ) {
							bestenergy() = currentenergy;
							best_state_on_node = state_on_node;
							if ( false ) { // hacking ///////////////////////////////////////////
								TR << "best-accept: ";
								for ( Size i=1; i<= Size(nmoltenres); ++i ) {
									if ( state_on_node( i ) == 0 ) {
										TR << '.';
									} else {
										RotamerSetCOP rotamer_set( rotamer_sets()->rotamer_set_for_moltenresidue( i ) );
										conformation::ResidueCOP rotamer( rotamer_set->rotamer( state_on_node( i ) ) );
										if ( rotamer->is_DNA() ) TR << rotamer->name1();
									}
								}
								TR << ' ' << nn << ' ' << n << ' ' << currentenergy << '\n';
							} // end hacking ///////////////////////////////////////
						}

					} else {
						// reject
						TR.Trace << "rejecting coupled rotamer substitution" << std::endl;
						ig_->consider_substitution( moltenres_id, prevrotamer_state, delta_energy, previous_energy_for_node );
						tmp_currentenergy = ig_->commit_considered_substitution();

						{ // debugging
							Real const dev( std::abs( tmp_currentenergy - currentenergy ) );
							if ( dev > 0.01 ) {
								TR << "equal3? " << dev << ' ' << tmp_currentenergy << ' ' << currentenergy << '\n';
							}
						} // scope
						currentenergy = tmp_currentenergy;
					} // accept or reject?

					loopenergy(nn) = currentenergy;

					debug_assert( !calc_rot_freq() );
					continue; // skip the logic below for single-rotamer substitution ////////////////////////////////////
				}
			}


			//std::cerr << "mres: " << moltenres_id << ", state: ";
			//std::cerr << rotamer_state_on_moltenres << ", deltaE: " << delta_energy << std::endl;

			//bk keep new rotamer if it is lower in energy or accept it at some
			//bk probability if it is higher in energy, if it is the first
			//bk rotamer to be tried at this position automatically accept it.
			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
				//std::cerr << "Accepted" << std::endl;
				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					bestenergy() = currentenergy;
					best_state_on_node = state_on_node;
				}
			} // accept

			loopenergy(nn) = currentenergy;
			core::PackerEnergy const temperature = get_temperature();

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

	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In FixbbCoupledRotamerSimAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
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
}

}//end namespace annealer
}//end namespace pack
}//end namespace core
