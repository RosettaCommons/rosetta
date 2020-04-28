// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/FixbbSimAnnealer.cc
/// @brief  This annealer (which is not being maintained or compiled) was used to generate the original training data for the smart annealer
/// @author Jack Maguire, jackmaguire1444@gmail.com

#define COLLECT_DATA

// Unit Headers
#include <core/pack/annealer/FixbbSimAnnealer.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack_basic/RotamerSetsBase.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <utility>
#include <utility/exit.hh>

#include <iostream>
#include <fstream>
#include <math.h>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using namespace ObjexxFCL;

static basic::Tracer TR( "core.pack.annealer.FixbbSimAnnealer" );

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
	ig_(std::move(ig)),
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
FixbbSimAnnealer::~FixbbSimAnnealer() = default;


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

	/// if pose has water molecules then check if virtualization is allowed
	/// (default=true) so that water molecules may be virtualized 50% of the time
	bool include_vrt = basic::options::option[ basic::options::OptionKeys::corrections::water::include_vrt ].value();

	std::ofstream annealer_trajectory;
	if ( record_annealer_trajectory_ ) {
		annealer_trajectory.open(trajectory_file_name_.c_str() );
	}

#ifdef COLLECT_DATA
	struct Key {
		core::Size mres;
		char AA;

		bool operator == ( Key const & o ) const {
			return o.mres == mres && o.AA == AA;
		}
	};

	//https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
	struct KeyHasher {
		std::size_t operator() ( Key const & k ) const {
			return ((std::hash< char >()( k.AA )
				^ (std::hash< Size >()( k.mres ) << 1)) >> 1);
		}
	};

	struct Value {
		core::Size num_accept_to_AA = 0;
		core::Size num_total_to_AA = 0;
		core::Real frac_to_AA() const {
			if ( num_total_to_AA == 0 ) return -1;
			return core::Real( num_accept_to_AA ) / core::Real ( num_total_to_AA );
		}

		core::Size num_accept_from_AA = 0;
		core::Size num_total_from_AA = 0;
		core::Real frac_from_AA() const {
			if ( num_total_from_AA == 0 ) return -1;
			return core::Real( num_accept_from_AA ) / core::Real ( num_total_from_AA );
		}

		core::Size num_accept_intra_AA = 0;
		core::Size num_total_intra_AA = 0;
		core::Real frac_intra_AA() const {
			if ( num_total_intra_AA == 0 ) return -1;
			return core::Real( num_accept_intra_AA ) / core::Real ( num_total_intra_AA );
		}

	};

	std::unordered_map< Key, utility::vector1< Value >, KeyHasher > my_data_map;
	my_data_map.max_load_factor( 0.1 );

	for ( core::Size rotid = 1; rotid <= rotamer_sets()->nrotamers(); ++rotid ) {
		Key key;
		key.mres = rotamer_sets()->moltenres_for_rotamer( rotid );
		key.AA = rotamer_sets()->rotamer( rotid )->name1();

		utility::vector1< Value > & vec = my_data_map[ key ];
		if ( vec.empty() ) {
			vec.resize( outeriterations );
		}
	}

	utility::vector1< core::Real > temps( outeriterations );
#endif

	//outer loop
	for ( int nn = 1; nn <= outeriterations; ++nn ) {
		setup_temperature(loopenergy,nn);
		temps[ nn ] = get_temperature();

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

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );

			/// removed const from rotamer_state_on_moltenres for code that virtualizes waters half the time
			int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// for waters, set to virtual 50% of the time
			if ( ( include_vrt ) && ( rotamer_sets()->rotamer_for_moltenres( moltenres_id, rotamer_state_on_moltenres )->type().is_virtualizable_by_packer() ) ) {  // don't want to do this for waters without a virtual state, like TP3
				core::Size const nrot = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
				/// last state is the virtual state
				if ( numeric::random::rg().uniform() < (static_cast<core::Real>(nrot)/2.0-1)/static_cast<core::Real>(nrot) ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
			}

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// initializing to zero but should be updated below.
			core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );

			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
				delta_energy, previous_energy_for_node);

#ifdef COLLECT_DATA
			Key k_to;
			Key k_from;
			if ( prevrotamer_state != 0 ) {
				k_to.mres = moltenres_id;
				k_to.AA = rotamer_sets()->rotamer_for_moltenres( moltenres_id, rotamer_state_on_moltenres )->name1();
				utility::vector1< Value > & inner_map_to = my_data_map[ k_to ];
				Value & val_to = inner_map_to[ nn ];

				k_from.mres = moltenres_id;
				k_from.AA = rotamer_sets()->rotamer_for_moltenres( moltenres_id, prevrotamer_state )->name1();
				utility::vector1< Value > & inner_map_from = my_data_map[ k_from ];
				Value & val_from = inner_map_from[ nn ];

				if ( k_from.AA == k_to.AA ) {
					++val_to.num_total_intra_AA;
				} else {
					++val_to.num_total_to_AA;
					++val_from.num_total_from_AA;
				}
			}

#endif

			//bk keep new rotamer if it is lower in energy or accept it at some
			//bk probability if it is higher in energy, if it is the first
			//bk rotamer to be tried at this position automatically accept it.
			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
#ifdef COLLECT_DATA
				if ( prevrotamer_state != 0 ) {
					utility::vector1< Value > & inner_map_to = my_data_map[ k_to ];
					Value & val_to = inner_map_to[ nn ];

					utility::vector1< Value > & inner_map_from = my_data_map[ k_from ];
					Value & val_from = inner_map_from[ nn ];

					if ( k_from.AA == k_to.AA ) {
						++val_to.num_accept_intra_AA;
					} else {
						++val_to.num_accept_to_AA;
						++val_from.num_accept_from_AA;
					}

				}
#endif

				//std::cout << " accepted\n";
				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					best_state_on_node = state_on_node;
					bestenergy() = currentenergy;
				}

				if ( record_annealer_trajectory_ ) {
					annealer_trajectory << moltenres_id << " " << rotamer_state_on_moltenres << " A\n";
				}

			} else if ( record_annealer_trajectory_ ) {
				annealer_trajectory << moltenres_id << " " << rotamer_state_on_moltenres << " R\n";
			}

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

#ifdef COLLECT_DATA
	//std::unordered_map< Key, utility::vector1< Value >, KeyHasher > my_data_map
	TR << "YYY_TEMP";
	for ( core::Real t : temps ) {
		TR << "," << t;
	}
	TR << std::endl;
	TR << "YYY_LNTEMP_plus0.1";
	for ( core::Real t : temps ) {
		TR << "," << log( t + 0.1 );
	}
	TR << std::endl;
	for ( std::pair< Key, utility::vector1< Value > > const & pair : my_data_map ) {
		Key const & k = pair.first;
		utility::vector1< Value > const & vals = pair.second;
		TR << "XXX_" << k.AA;
		core::Size count = 1;
		for ( Value const & v : vals ) {
			TR << "," << log( temps[ count ] + 0.1 ) << "," << v.frac_to_AA() << "," << v.frac_from_AA() << "," << v.frac_intra_AA();
			++count;
		}
		runtime_assert( int(k.mres) <= nmoltenres );
		auto const rotset = rotamer_sets()->rotamer_set_for_moltenresidue( k.mres );
		runtime_assert( rotset != nullptr );
		runtime_assert( best_state_on_node( k.mres ) != 0 );
		auto const rotamer = rotset->rotamer( best_state_on_node( k.mres ) );
		runtime_assert( rotamer != nullptr );
		char const final_AA = rotamer->name1();
		if ( final_AA == k.AA ) {
			TR << ",1.0";
		} else {
			TR << ",0.0";
		}
		TR << std::endl;
	}
#endif

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( int ii = 1; ii <= nmoltenres; ++ii ) {
		int const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}
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
