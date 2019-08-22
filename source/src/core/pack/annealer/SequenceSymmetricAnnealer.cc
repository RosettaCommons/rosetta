// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/SequenceSymmetricAnnealer.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com


#include <basic/Tracer.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/annealer/SequenceSymmetricAnnealer.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
//#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/cxx_versioning_macros.hh>
//C++
#include <fstream>
#include <iostream>
#include <utility>
#include <unordered_set>


using namespace ObjexxFCL;

//#ifndef NDEBUG
static basic::Tracer TR( "core.pack.annealer.SequenceSymmetricAnnealer" );
//#endif

namespace core {
namespace pack {
namespace annealer {

namespace {

core::Size
determine_num_chains( core::pose::Pose const & pose ){
	//Okay this implementation is kinda lazy, but is cold so it's okay
	//we need to count the number of chains that have nonvirtual residues in them
	//The goal is to avoid caring about virtual roots and similar ideas
	std::unordered_set< core::Size > chains;
	for ( core::Size resid = 1; resid <= pose.size(); ++resid ) {
		if ( ! pose.residue( resid ).is_virtual_residue() ) {
			chains.insert( pose.chain( resid ) );
		}
	}
	//TR << "chains.size(): " << chains.size() << std::endl;
	return chains.size();
}

}

////////////////////////////////////////////////////////////////////////////////
SequenceSymmetricAnnealer::SequenceSymmetricAnnealer(
	core::pose::Pose const & pose,
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
	starting_sequence_( pose.sequence() ),
	pdb_info_( pose.pdb_info() ),
	ig_(std::move(ig)),
	record_annealer_trajectory_( false )
{
	num_chains_ = determine_num_chains( pose );
}

SequenceSymmetricAnnealer::SequenceSymmetricAnnealer(
	core::pose::Pose const & pose,
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
	starting_sequence_( pose.sequence() ),
	pdb_info_( pose.pdb_info() ),
	ig_(std::move(ig)),
	record_annealer_trajectory_( false )
{
	num_chains_ = determine_num_chains( pose );
}

/// @brief virtual destructor
SequenceSymmetricAnnealer::~SequenceSymmetricAnnealer() = default;

namespace {

//TODO measure benchmarks using utility::vector1< std::basic_string< Size > > for small string optimizations
NODISCARD
utility::vector1< utility::vector1< Size > >
create_corresponding_mress_for_mres(
	Size const nmoltenres,
	pose::PDBInfoCOP const pdb_info,
	rotamer_set::FixbbRotamerSets const & rotamer_sets
){
	utility::vector1< utility::vector1< Size > > corresponding_mress_for_mres( nmoltenres );

	using Num  = Size;
	using Mres = Size;

	std::unordered_map< Num, std::list< Mres > > mress_for_num;
	mress_for_num.max_load_factor( 0.1 );

	for ( Size mres = 1; mres <= nmoltenres; ++mres ) {
		Size const resid = rotamer_sets.moltenres_2_resid( mres );
		Size const pdb_num = pdb_info->number( resid );

		mress_for_num[ pdb_num ].push_back( mres );
	}

#ifndef NDEBUG
	TR << "START RESID CLUSTERS" << std::endl;
#endif

	for ( std::pair< Num, std::list< Mres > > const & iter : mress_for_num ) {
		std::list< Mres > const & mress = iter.second;

		for ( auto iter = mress.begin(); iter != mress.end(); ++iter ) {
			Mres const mres1 = * iter;
#ifndef NDEBUG
			TR << mres1 << " ";
#endif
			for ( auto iter2 = std::next( iter ); iter2 != mress.end(); ++iter2 ) {
				Mres const mres2 = * iter2;
				corresponding_mress_for_mres[ mres1 ].push_back( mres2 );
				corresponding_mress_for_mres[ mres2 ].push_back( mres1 );
				//TR << "ADDING " << mres1 << " and " << mres2 << std::endl;
			}
		}
#ifndef NDEBUG
		TR << std::endl;
#endif

	}

#ifndef NDEBUG
	TR << "END RESID CLUSTERS" << std::endl;
#endif


	return corresponding_mress_for_mres;
}

}

/// @brief sim_annealing for fixed backbone design mode
void SequenceSymmetricAnnealer::run()
{
	Size const nmoltenres = ig_->get_num_nodes();

	utility::vector1< utility::vector1< Size > > corresponding_mress_for_mres =
		create_corresponding_mress_for_mres( nmoltenres, pdb_info_, * rotamer_sets() );

#ifndef NDEBUG
	TR << "Linked Resids:" << std::endl;
	for ( core::Size mres = 1; mres <= corresponding_mress_for_mres.size(); ++mres ) {
		TR << "mres_" << mres;
		auto const & v = corresponding_mress_for_mres[ mres ];
		for ( auto const i : v ) {
			TR << " " << i;
		}
		TR << std::endl;
	}
	TR << "End Linked Resids" << std::endl;
#endif

	FArray1D_int state_on_node( nmoltenres,0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D_float loopenergy(maxouteriterations,0.0);

	//bk variables for calculating rotamer frequencies during simulation
	//int nsteps = 0;
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

	/// if pose has water molecules then check if virtualization is allowed
	/// (default=true) so that water molecules may be virtualized 50% of the time
	bool include_vrt = basic::options::option[ basic::options::OptionKeys::corrections::water::include_vrt ].value();

	std::ofstream annealer_trajectory;
	if ( record_annealer_trajectory_ ) {
		annealer_trajectory.open(trajectory_file_name_.c_str() );
	}

	//outer loop
	for ( int nn = 1; nn <= outeriterations; ++nn ) {
		setup_temperature(loopenergy,nn);
		if ( quench() ) {
			currentenergy = bestenergy();
			state_on_node = best_state_on_node;
			ig_->set_network_state( state_on_node );
		}

		int inneriterations = get_inneriterations();

		float treshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( treshold_for_deltaE_inaccuracy );

		//inner loop
		for ( int n = 1; n <= inneriterations; ++n ) {
			int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			char const name1 = rotamer_sets()->rotamer( ranrotamer )->name1();
			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int const resid = rotamer_sets()->moltenres_2_resid( moltenres_id );

			//Do not let position mutate if we can not mutate its partner
			if ( corresponding_mress_for_mres[ moltenres_id ].size() != num_chains_ - 1 ) {
				//then there is at least one symmetric partner that cannot mutate/repack

				//TR << "corresponding_mress_for_mres[ moltenres_id ].size(): " << corresponding_mress_for_mres[ moltenres_id ].size() << std::endl;

				//only allow if this does not cause a mutation
				if ( name1 != starting_sequence_[ resid - 1 ] ) {
					continue;
				}
			}


			/// removed const from rotamer_state_on_moltenres for code that virtualizes waters half the time
			int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// for waters, set to virtual 50% of the time
			core::conformation::Residue curres( *rotamer_sets()->rotamer_for_moltenres(moltenres_id, rotamer_state_on_moltenres) );
			if ( ( curres.name3() == "HOH" ) && ( include_vrt ) ) {  // don't want to do this for waters without a virtual state, like TP3
				auto const rand_num = numeric::random::rg().uniform();
				auto const nrot = rotamer_sets()->nrotamers_for_moltenres( moltenres_id );
				/// last state is the virtual state
				if ( rand_num < (nrot/2.0-1)/nrot ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
			}

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			if ( corresponding_mress_for_mres.size() == 0 ) {
				//Do normal protocol

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

				continue;

			} else {
				//Do fancy protocol
				struct SavedState {
					SavedState() = default;
					SavedState( int m, int p ) :
						mres( m ),
						previous_state( p )
					{}

					int mres;
					int previous_state;
				};

				bool any_previous_unassigned = false;

				std::list< SavedState > starting_state;
				starting_state.emplace_back( moltenres_id, prevrotamer_state );
				any_previous_unassigned |= ( prevrotamer_state == 0 );
				bool fail = false;
				for ( auto const other_moltenres_id : corresponding_mress_for_mres[ moltenres_id ] ) {
					//cache current state
					auto const state = state_on_node( other_moltenres_id );
					any_previous_unassigned |= ( state == 0 );
					starting_state.emplace_back( other_moltenres_id, state );

					//ensure that there is at least one candidate state
					rotamer_set::RotamerSetCOP other_rotamer_set =
						rotamer_sets()->rotamer_set_for_moltenresidue( other_moltenres_id );
					Size const num_other_rots = other_rotamer_set->num_rotamers();
					bool at_least_one_good_one = false;
					//Size num_good_ones = 0;
					for ( Size other_rot_id=1; other_rot_id<=num_other_rots; ++other_rot_id ) {
						if ( other_rotamer_set->rotamer( other_rot_id )->name1() == name1 ) {
							at_least_one_good_one = true;
							break;
							//++num_good_ones;
						}
					}
					fail |= ! at_least_one_good_one;
					//TR << num_good_ones << " good rots for " << other_moltenres_id << std::endl;
				}

				if ( fail ) continue;

				core::PackerEnergy global_previous_energy = 0;
				core::PackerEnergy global_deltaE = 0;

				core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );
				ig_->consider_substitution(
					//input:
					moltenres_id,
					rotamer_state_on_moltenres,
					//output:
					delta_energy,
					previous_energy_for_node
				);
				auto current_energy = ig_->commit_considered_substitution();

				global_previous_energy = previous_energy_for_node;
				global_deltaE = delta_energy;

				std::map< int /*Node*/, int /*state*/ > state_for_node_map;

				for ( auto const other_moltenres_id : corresponding_mress_for_mres[ moltenres_id ] ) {

					//Determine candidate rotamers for this position
					rotamer_set::RotamerSetCOP other_rotamer_set =
						rotamer_sets()->rotamer_set_for_moltenresidue( other_moltenres_id );
					Size const num_other_rots = other_rotamer_set->num_rotamers();
					utility::vector1< Size > local_ids_for_good_rotamers;
					for ( Size other_rot_id=1; other_rot_id<=num_other_rots; ++other_rot_id ) {
						if ( other_rotamer_set->rotamer( other_rot_id )->name1() == name1 ) {
							local_ids_for_good_rotamers.push_back( other_rot_id );
						}
					}

					//TR << local_ids_for_good_rotamers.size() << " good rots for " << other_moltenres_id << std::endl;


					//This should have been checked ~30 lines ago
					runtime_assert( ! local_ids_for_good_rotamers.empty() );//failure here is failure of the code, not the runtime

					//Pick a good rotamer
					auto const rand_index =
						numeric::random::random_range( 1, local_ids_for_good_rotamers.size() );
					int const other_rotamer_state = local_ids_for_good_rotamers[ rand_index ];
					ig_->consider_substitution(
						//input:
						other_moltenres_id,
						other_rotamer_state,
						//output:
						delta_energy,
						previous_energy_for_node
					);
					current_energy = ig_->commit_considered_substitution();
					global_previous_energy += previous_energy_for_node;
					global_deltaE += delta_energy;

					state_for_node_map[ other_moltenres_id ] = other_rotamer_state;
				}

				core::Size const num_changes =
					corresponding_mress_for_mres[ moltenres_id ].size() + 1;
				core::PackerEnergy const previous_energy_average =
					global_previous_energy / num_changes;
				core::PackerEnergy const delta_energy_average = global_deltaE / num_changes;

				if ( any_previous_unassigned ||
						pass_metropolis( previous_energy_average, delta_energy_average ) ) {

					state_on_node( moltenres_id ) = rotamer_state_on_moltenres;
					for ( auto const & node_state_pair : state_for_node_map ) {
						state_on_node( node_state_pair.first ) = node_state_pair.second;
					}

					currentenergy = current_energy;
					if ( any_previous_unassigned || (currentenergy < bestenergy() ) ) {
						best_state_on_node = state_on_node;
						bestenergy() = currentenergy;
					}

				} else { //reject
					for ( SavedState const state : starting_state ) {
						ig_->consider_substitution(
							//input:
							state.mres,
							state.previous_state,
							//output:
							delta_energy,
							previous_energy_for_node
						);
						currentenergy = ig_->commit_considered_substitution();
					}
				}

			}

			loopenergy(nn) = currentenergy;
			//float const temperature = get_temperature();

		} // end of inneriteration loop
	} //end of outeriteration loop

	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In SequenceSymmetricAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
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
	for ( uint ii = 1; ii <= nmoltenres; ++ii ) {
		int const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}
}

void SequenceSymmetricAnnealer::record_annealer_trajectory( bool setting ) {
	record_annealer_trajectory_ = setting;
}

void SequenceSymmetricAnnealer::trajectory_file_name( std::string const & setting ) {
	trajectory_file_name_ = setting;
}


}//end namespace annealer
}//end namespace pack
}//end namespace core
