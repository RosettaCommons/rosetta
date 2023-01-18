// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/SequenceSymmetricAnnealer.cc
/// @author Jack Maguire
/// @author Updated by Tim Neary, timdot10@gmail.com


#include <basic/Tracer.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/annealer/SequenceSymmetricAnnealer.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/select/residue_selector/CachedResidueSubset.hh>

#include <numeric/random/random.hh>

#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// External headers

//C++
#include <fstream>
#include <iostream>
#include <utility>
#include <unordered_set>

#include <ObjexxFCL/FArray1D.hh> // AUTO IWYU For FArray1D


using namespace ObjexxFCL;

//#ifndef NDEBUG
static basic::Tracer TR( "core.pack.annealer.SequenceSymmetricAnnealer" );
//#endif

namespace core {
namespace pack {
namespace annealer {


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
	FArray1D_float & rot_freq,
	std::string const & rs_prefix
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
	record_annealer_trajectory_( false ),
	rs_prefix_( rs_prefix )
{
	search_pose_for_residue_subsets( pose );
	power_mode_ = ! residue_subsets_.empty();
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
	FArray1D_float & rot_freq,
	std::string const & rs_prefix
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
	record_annealer_trajectory_( false ),
	rs_prefix_( rs_prefix )
{
	search_pose_for_residue_subsets( pose );
	power_mode_ = ! residue_subsets_.empty();
}

/// @brief virtual destructor
SequenceSymmetricAnnealer::~SequenceSymmetricAnnealer() = default;

core::Size
process_linked_res_set(
	std::unordered_set< core::Size> & set,
	utility::vector1< utility::vector1< Size > > const & corresponding_mress_for_mres,
	utility::vector1< bool > & processed ) {
	std::unordered_set< core::Size > set_copy(set);
	for ( core::Size const lkd_mres : set_copy ) { // now create an union from all linked mres
		if ( processed[ lkd_mres ] == true ) {
			continue; // skip if we have already processed this mres.
		}

		for ( core::Size jj = 1; jj <=corresponding_mress_for_mres[ lkd_mres ].size(); ++jj ) {
			set.emplace( corresponding_mress_for_mres[ lkd_mres ][ jj ] );
		}
		processed[ lkd_mres ] = true;
	}
	return set_copy.size();
}

//NODISCARD_ATTR
utility::vector1< utility::vector1< Size > >
SequenceSymmetricAnnealer::create_corresponding_mress_for_mres() const {
	struct Key {
		core::Size idx;
		core::Size sele_uid = 0;

		bool operator == ( Key const & o ) const {
			return o.idx == idx && o.sele_uid == sele_uid;
		}
	};

	//https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
	struct KeyHasher {
		std::size_t operator() ( Key const & k ) const {
			return ((std::hash< core::Size >()( k.idx )
				^ (std::hash< core::Size >()( k.sele_uid ) << 1)) >> 1);
		}
	};

	using SizePair = std::pair< core::Size, core::Size >;
	struct SizePairHash { // So pair can be hased
		std::size_t operator() ( std::pair<core::Size, core::Size> const & pair ) const
		{
			return (std::hash< core::Size >()( pair.first )
				^ std::hash< core::Size >()( pair.second ));
		}
	};

	Size const nmoltenres = ig_->get_num_nodes();
	utility::vector1< utility::vector1< Size > > corresponding_mress_for_mres( nmoltenres );

	std::unordered_map< Key, std::list< core::Size >, KeyHasher > idx_for_key;
	idx_for_key.max_load_factor( 0.1 );
	std::unordered_map< SizePair, core::Size, SizePairHash > idx_tracker;
	idx_tracker.max_load_factor( 0.1 );

	if ( residue_subsets_.size() == 0 ) return corresponding_mress_for_mres; // We have no links so want to return empty vect of vect
	core::Size num_res = residue_subsets_[ 1 ][ 1 ].size(); // Technically this may fail if the corresponding residue selector is cached and then residues are added/removed

	for ( core::Size resid = 1; resid <= num_res; ++resid ) {
		Key key;
		for ( core::Size linked_subset = 1; linked_subset <= residue_subsets_.size(); ++linked_subset ) {
			key.sele_uid = linked_subset;

			for ( core::Size sele_id_in_group = 1; sele_id_in_group <= residue_subsets_[ linked_subset ].size(); ++sele_id_in_group ) {
				if ( residue_subsets_[ linked_subset ][ sele_id_in_group ][ resid ] ) {
					// std::unordered_map< SizePair, core::Size >::iterator
					auto it = idx_tracker.find( SizePair( linked_subset, sele_id_in_group ) ); // get iterator to key
					core::Size idx;
					if ( it == idx_tracker.end() ) {
						idx = 0;
						idx_tracker[ SizePair( linked_subset, sele_id_in_group ) ] = 0;
					} else {
						it->second += 1;
						idx = it->second;
					}
					key.idx = idx;
					idx_for_key[ key ].emplace_back( resid );
				}
			}
		}
	}

	std::unordered_map< core::Size, core::Size > mres_to_resid; // Using map as pdb_num may be greater than nmoltenres

	// get map of resid -> mres
	for ( Size mres = 1; mres <= nmoltenres; ++mres ) {
		//TR << mres << " " << rotamer_sets()->nrotamers_for_moltenres( mres ) << std::endl;

		if ( rotamer_sets()->nrotamers_for_moltenres( mres ) == 0 ) continue;
		if ( rotamer_sets()->rotamer_for_moltenres( mres, 1 )->is_virtual_residue() ) continue;
		Size const resid = rotamer_sets()->moltenres_2_resid( mres );
		// Size const pdb_num = pdb_info_->number( resid );

		mres_to_resid[ resid ] = mres;
	}

	// Build linked residues using mres not resids
	for ( std::pair< const Key, std::list< core::Size > > const & iter : idx_for_key ) {
		std::list< core::Size > const & res_ids = iter.second;

		for ( auto iter1 = res_ids.begin(); iter1 != res_ids.end(); ++iter1 ) {
			auto mres1 = mres_to_resid.find( *iter1 );
			if (  mres1 == mres_to_resid.end() ) continue;

			for ( auto iter2 = std::next( iter1 ); iter2 != res_ids.end(); ++iter2 ) {
				auto mres2 = mres_to_resid.find( *iter2 );
				if (  mres2 == mres_to_resid.end() ) {
					continue;
				}
				corresponding_mress_for_mres[ mres1->second ].push_back( mres2->second );
				corresponding_mress_for_mres[ mres2->second ].push_back( mres1->second );
				//TR << "ADDING " << mres1 << " and " << mres2 << std::endl;
			}
		}
	}

	// Now ensure that we have union of all linked residues and remove duplicate entries
	utility::vector1< bool > processed( corresponding_mress_for_mres.size(), false ); // used to speed up processing
	for ( core::Size ii = 1; ii <= corresponding_mress_for_mres.size(); ++ii ) {
		if ( processed[ ii ] == true ) continue; // Skip this mres if we have already processed it

		std::unordered_set< core::Size > set;
		set.max_load_factor( 0.1 );
		for ( core::Size jj = 1; jj <= corresponding_mress_for_mres[ ii ].size(); ++jj ) {
			set.emplace( corresponding_mress_for_mres[ ii ][ jj ] );
		}
		processed[ ii ] = true;

		core::Size count = 0;
		while ( count < set.size() ) {
			// need to repeat the set creation until all res are linked
			// This process can almost certainly be optimised but is only ever done once so its probably ok as is.
			count = process_linked_res_set( set, corresponding_mress_for_mres, processed );
		}

		for ( core::Size const mres : set ) {
			utility::vector1< core::Size > & vec = corresponding_mress_for_mres[ mres ];
			vec.clear(); // Need to clear current contents as we do not want duplicates
			vec.resize( set.size() - 1 );
			std::copy_if( set.begin(), set.end(), vec.begin(),
				[=]( core::Size i ){ return i != mres; } );
			// copy everything but the mres itself
			processed[ mres ] = true; // Should help to speed up creation of linked mres as all linked res will also be removed from the list
		}
	}
	return corresponding_mress_for_mres;
}

void
SequenceSymmetricAnnealer::update_shared_residue_map(
	core::Size const id,
	std::unordered_map< std::string, core::Size > & map ) const {

	rotamer_set::RotamerSetCOP rotamer_set =
		rotamer_sets()->rotamer_set_for_moltenresidue( id );
	Size const num_rots = rotamer_set->num_rotamers();

#ifndef NDEBUG
	TR << "Residue type set for resid (Rosetta numbering): " << std::to_string( rotamer_sets()->moltenres_2_resid( id ) ) << " = ";
#endif
	std::unordered_set< std::string > curr_set; // To store all res_types
	curr_set.max_load_factor( 0.1 );
	// Get the current set of residue types for the given mres (id)
	for ( Size rot_id = 1; rot_id <= num_rots; ++rot_id ) {
		std::string const & res_type = rotamer_set->rotamer( rot_id )->type().interchangeability_group();
		curr_set.insert( res_type );
	}

	// Map contains a count of the number of times a residue type was seen for each linked mres
	// now we need to update the map with whether the res type was seen for the current mres
	for ( std::string const & res_n1 : curr_set ) {
#ifndef NDEBUG
		TR << res_n1 << ", ";
#endif
		if ( map.find( res_n1 ) == map.end() ) {
			map.emplace( std::make_pair( res_n1, 1 ) );
		} else {
			map[ res_n1 ] += 1;
		}
	}
#ifndef NDEBUG
	TR << std::endl;
#endif
}

std::unordered_set< std::string >
SequenceSymmetricAnnealer::get_shared_residue_types(
	core::Size const curr_mres,
	utility::vector1< core::Size > const & linked_res ) const {

	std::unordered_set< std::string > common_res_types;
	common_res_types.max_load_factor( 0.1 );

	if ( linked_res.size() == 0 ) {
		return common_res_types;
	}
	// A map with a count of the number of times a residue type was seen for a set of mres
	std::unordered_map< std::string, core::Size > shared_res_types;

	// Now we need to get a set of the residue types present for a given set of linked_res
	// As linked_res does not include the current res, must add it here.
	update_shared_residue_map( curr_mres, shared_res_types );

	for ( core::Size ii = 1; ii <= linked_res.size(); ++ii ) {
		update_shared_residue_map( linked_res[ ii ], shared_res_types ); // Update map for the linked mres

		// Determine if any (rotamer) residue types are shared across all linked res
		bool any_shared = false;
		for ( auto const & it : shared_res_types ) {
			if ( it.second == ii + 1 ) { // Accounts for addition of current res too
				any_shared = true;
				break;
			}
		}
		if ( ! any_shared ) {
			return common_res_types;
		}
	}

	// Now make set with all res types found in each linked residue.
	for ( auto const & it : shared_res_types ) {
		if ( it.second == linked_res.size() + 1 ) common_res_types.insert( it.first );
	}
	return common_res_types;
}

void
SequenceSymmetricAnnealer::print_linked_residues( utility::vector1< utility::vector1< core::Size > > const & corresponding_mress_for_mres ) const {
	TR << "Linked Resids (Rosetta numbering):" << std::endl;
	for ( core::Size mres = 1; mres <= corresponding_mress_for_mres.size(); ++mres ) {
		Size const resid_1 = rotamer_sets()->moltenres_2_resid( mres );
		TR << "Residue " << resid_1 << " links: ";
		auto const & v = corresponding_mress_for_mres[ mres ];
		for ( auto const i : v ) {
			Size const resid_2 = rotamer_sets()->moltenres_2_resid( i );
			TR << " " << resid_2;
		}
		TR << std::endl;
	}
	TR << "End Linked Resids" << std::endl;
}

SeqSymmAnnealerSetup
SequenceSymmetricAnnealer::setup_for_linked_residues() {
	utility::vector1< utility::vector1< Size > > const corresponding_mress_for_mres =
		create_corresponding_mress_for_mres();
	print_linked_residues( corresponding_mress_for_mres );

	// Determine if any mres do not share any residue types as this is an impossible problem for this annealer.
	utility::vector1< StringSetOP > common_res_types( ig_->get_num_nodes() );
	for ( core::Size ii = 1; ii <= corresponding_mress_for_mres.size(); ++ii ) {
		auto new_set = utility::pointer::make_shared< std::unordered_set< std::string > >(
			get_shared_residue_types( ii, corresponding_mress_for_mres[ ii ] ) );
		if ( ! corresponding_mress_for_mres[ ii ].empty() && new_set->empty() ) {
			// If there are linked res but no common residue types
			std::string exit_msg;
			exit_msg +=
				"Some linked residues were found to not share any residue types, though there may be more. "
				"The linked residues (Rosetta numbering) were: " + std::to_string( rotamer_sets()->moltenres_2_resid( ii ) ) + " ";
			for ( core::Size const mres : corresponding_mress_for_mres[ ii ] ) exit_msg += std::to_string( rotamer_sets()->moltenres_2_resid( mres ) ) + " ";
			utility_exit_with_message( exit_msg );
		}

		for ( auto const & mres : corresponding_mress_for_mres[ ii ] ) common_res_types[ mres ] = new_set; // assign all linked res the common set.
	}

	return SeqSymmAnnealerSetup{ corresponding_mress_for_mres, common_res_types };
}

/// @brief sim_annealing for fixed backbone design mode
void
SequenceSymmetricAnnealer::run() {
	auto setup_info = setup_for_linked_residues(); // Setup for annealer run, get linked res information.
	utility::vector1< utility::vector1< core::Size > > corresponding_mress_for_mres =
		std::move( setup_info.corresponding_mress_for_mres );
	utility::vector1< StringSetOP > common_res_types = // StringSetOP == utility::pointer::shared_ptr< std::unordered_set< std::string > >
		std::move( setup_info.common_res_types );

	Size const nmoltenres = ig_->get_num_nodes();

	FArray1D_int state_on_node( nmoltenres,0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D_float loopenergy(maxouteriterations,0.0);

	//bk variables for calculating rotamer frequencies during simulation
	//int nsteps = 0;
	FArray1D_int nsteps_for_rot( ig_->get_num_total_states(), 0 );

	//--------------------------------------------------------------------
	//initialize variables

	core::PackerEnergy currentenergy = 0.0;

	ig_->prepare_graph_for_simulated_annealing();
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

		float threshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( threshold_for_deltaE_inaccuracy );

		//inner loop
		for ( int n = 1; n <= inneriterations; ++n ) {
			int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			std::string const & interchangeability_group = rotamer_sets()->rotamer( ranrotamer )->type().interchangeability_group();
			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			// int const resid = rotamer_sets()->moltenres_2_resid( moltenres_id );

			/// removed const from rotamer_state_on_moltenres for code that virtualizes waters half the time
			int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state && corresponding_mress_for_mres[ moltenres_id ].empty() ) continue;
			// can only implicitly know to skip iteration if there is no links

			// for waters, set to virtual 50% of the time
			core::conformation::Residue curres( *rotamer_sets()->rotamer_for_moltenres(moltenres_id, rotamer_state_on_moltenres) );
			if ( ( curres.name3() == "HOH" ) && ( include_vrt ) ) {  // don't want to do this for waters without a virtual state, like TP3
				auto const rand_num = numeric::random::rg().uniform();
				auto const nrot = rotamer_sets()->nrotamers_for_moltenres( moltenres_id );
				/// last state is the virtual state
				if ( rand_num < (nrot/2.0-1)/nrot ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
			}

			if ( corresponding_mress_for_mres[ moltenres_id ].empty() ) {
				//Do normal protocol

				core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );
				ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
					delta_energy, previous_energy_for_node);

				if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
					currentenergy = ig_->commit_considered_substitution();
					state_on_node(moltenres_id) = rotamer_state_on_moltenres;
					if ( (prevrotamer_state == 0)||( currentenergy < bestenergy() ) ) {
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

				// If picked a non common residue type then skip iteration
				if ( common_res_types[ moltenres_id ]->find( interchangeability_group ) == common_res_types[ moltenres_id ]->end() ) continue;

				bool any_previous_unassigned = false;
				utility::vector1< SavedState > starting_state;

				starting_state.emplace_back( moltenres_id, prevrotamer_state );
				any_previous_unassigned |= ( prevrotamer_state == 0 );

				for ( auto const other_moltenres_id : corresponding_mress_for_mres[ moltenres_id ] ) {
					//cache current state
					auto const state = state_on_node( other_moltenres_id );
					any_previous_unassigned |= ( state == 0 );
					starting_state.emplace_back( other_moltenres_id, state );
				}

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
				core::Size num_changes = 0;

				for ( core::Size ii = 1; ii <= corresponding_mress_for_mres[ moltenres_id ].size(); ++ii ) {
					core::Size const other_moltenres_id = corresponding_mress_for_mres[ moltenres_id ][ ii ];
					//Determine candidate rotamers for this position
					rotamer_set::RotamerSetCOP other_rotamer_set =
						rotamer_sets()->rotamer_set_for_moltenresidue( other_moltenres_id );
					Size const num_other_rots = other_rotamer_set->num_rotamers();
					utility::vector1< Size > local_ids_for_good_rotamers;
					for ( Size other_rot_id = 1; other_rot_id <= num_other_rots; ++other_rot_id ) {
						if ( other_rotamer_set->rotamer( other_rot_id )->type().interchangeability_group() == interchangeability_group ) {
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

					if ( starting_state[ ii ].previous_state != other_rotamer_state ) {
						// Can skip the subsitution consideriation if the new random state is identical to previous.

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
						++num_changes;
					}
					state_for_node_map[ other_moltenres_id ] = other_rotamer_state;
				}

				core::PackerEnergy const previous_energy_average =
					global_previous_energy / num_changes;
				core::PackerEnergy const delta_energy_average = global_deltaE / num_changes;

				if ( any_previous_unassigned ||
						pass_metropolis( previous_energy_average, delta_energy_average ) ) {
					// Accept changes if any state is previously unassigned or the average change in energy passes the MC.

					state_on_node( moltenres_id ) = rotamer_state_on_moltenres;
					for ( auto const & node_state_pair : state_for_node_map ) {
						state_on_node( node_state_pair.first ) = node_state_pair.second;
					}

					currentenergy = current_energy;
					if ( any_previous_unassigned || ( currentenergy < bestenergy() ) ) {
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

void
SequenceSymmetricAnnealer::record_annealer_trajectory( bool setting ) {
	record_annealer_trajectory_ = setting;
}

void
SequenceSymmetricAnnealer::trajectory_file_name( std::string const & setting ) {
	trajectory_file_name_ = setting;
}

//MUCH of this is infuenced by StoredResidueSubsetSelector.cc
void
SequenceSymmetricAnnealer::search_pose_for_residue_subsets( core::pose::Pose const & pose ) {
	residue_subsets_.clear();

	//Look for selections
	if ( ! pose.data().has( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET ) ) {
		return;
	}

	auto const temp_ptr = pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET );
	core::select::residue_selector::CachedResidueSubset const & stored_subsets =
		*( utility::pointer::static_pointer_cast< core::select::residue_selector::CachedResidueSubset const >( temp_ptr ) );

	for ( core::Size ii = 0; true; ++ii ) {
		if ( ! stored_subsets.has_subset_prefix( rs_prefix_ + std::to_string(ii) ) ) {
			break;
		}
		utility::vector1< utility::vector1< bool > > linked_subsets;
		for ( core::Size jj = 0; true; ++jj ) {
			std::string const magic_selector_name =
				rs_prefix_ + std::to_string(ii) + "_" + std::to_string(jj);
			if ( ! stored_subsets.has_subset( magic_selector_name ) ) {
#ifndef NDEBUG
				TR << "Subset " << magic_selector_name << " is absent" << std::endl;
#endif
				break;
			}
			core::select::residue_selector::ResidueSubsetCOP subset =
				stored_subsets.get_subset( magic_selector_name );
			runtime_assert( subset != nullptr );
#ifndef NDEBUG
			TR << "Printing Subset # " << ii << "," << jj << std::endl;
			for ( core::Size resid = 1; resid <= subset->size(); ++resid ) {
				TR << resid << " " << (*subset)[ resid ] << std::endl;
			}
#endif
			linked_subsets.emplace_back( *subset );
		}
		residue_subsets_.emplace_back( linked_subsets );
	}

}

}//end namespace annealer
}//end namespace pack
}//end namespace core
