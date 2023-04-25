// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/lddt.cc
/// @brief Utilities for calcuating lDDT.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/scoring/lddt.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/automorphism.hh>
#include <core/id/AtomID.hh>

#include <utility/LexicographicalIterator.hh>
#include <basic/Tracer.hh>

#include <boost/container/flat_map.hpp>

namespace core {
namespace scoring {

static basic::Tracer TR("core.scoring.lddt");


///////////////////////////////////////////////

core::Real
lddt(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	bool consider_alt
) {
	lDDT_Calculator calc;
	calc.consider_alt_states( consider_alt );

	return calc( ref, model );
}

utility::vector1< core::Real >
per_res_lddt(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	bool consider_alt
) {
	lDDT_Calculator calc;
	calc.consider_alt_states( consider_alt );

	return calc.residue_lDDT( ref, model );
}

//////////// lDDT_Calculator ///////////////////////////////

lDDT_Calculator::lDDT_Calculator( PredicateCOP pred ) {
	predicate( pred );
}

void
lDDT_Calculator::predicate( PredicateCOP pred ) {
	if ( pred ) {
		predicate_ = pred;
	} else {
		predicate_ = utility::pointer::make_shared< IsHeavyAtomPredicate  >();
	}
}


core::Real
lDDT_Calculator::operator() (
	core::pose::Pose const & ref,
	core::pose::Pose const & model
) const {
	runtime_assert( ref.size() == model.size() );
	return operator()(ref, model, make_identity_map(ref) );
}

core::Real
lDDT_Calculator::operator() (
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Size > const & res_map
) const {
	std::map< core::Size, core::Real > res;
	std::map< core::id::AtomID, core::Real > atom;

	return get_stats(ref, model, res_map, false, res, false, atom );
}

utility::vector1< core::Real >
lDDT_Calculator::residue_lDDT(
	core::pose::Pose const & ref,
	core::pose::Pose const & model
) const {
	runtime_assert( ref.size() == model.size() );
	std::map< core::Size, core::Real > res;

	residue_lDDT(ref, model, make_identity_map(ref), res);

	utility::vector1< core::Real > residue_lddt;
	for ( core::Size ii(1); ii <= ref.size(); ++ii ) {
		residue_lddt.push_back( res[ii] ); // Defaults to zero if not present in map.
	}

	return residue_lddt;
}

void
lDDT_Calculator::residue_lDDT(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Size > const & res_map,
	std::map< core::Size, core::Real > & residue_lddt // return by reference
) const {
	std::map< core::id::AtomID, core::Real > atom;

	get_stats( ref, model, res_map, true, residue_lddt, false, atom );
}

/// @brief Get the atom lDDTs for all atoms in the poses
std::map< core::id::AtomID, core::Real >
lDDT_Calculator::atom_lDDT(
	core::pose::Pose const & ref,
	core::pose::Pose const & model
) const {
	runtime_assert( ref.size() == model.size() );

	std::map< core::id::AtomID, core::Real > atom_lddt;
	atom_lDDT( ref, model, make_identity_map(ref), atom_lddt );
	return atom_lddt;
}

/// @brief Get the atom lDDTs, for the residue pairings specified in res_map.
/// Only reference residues in res_map are considered for pairings.
/// To include residues but set them unmatched, map them to zero.
void
lDDT_Calculator::atom_lDDT(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Size > const & res_map,
	std::map< core::id::AtomID, core::Real > & atom_lddt // return by reference
) const {
	std::map< core::Size, core::Real > res;

	get_stats( ref, model, res_map, false, res, true, atom_lddt );
}

core::Real
lDDT_Calculator::all_lDDT(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Real > & residue_lddt, // return by reference
	std::map< core::id::AtomID, core::Real > & atom_lddt // return by reference
) const {
	runtime_assert( ref.size() == model.size() );
	return all_lDDT( ref, model, make_identity_map(ref), residue_lddt, atom_lddt );
}

core::Real
lDDT_Calculator::all_lDDT(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Size > const & res_map,
	std::map< core::Size, core::Real > & residue_lddt, // return by reference
	std::map< core::id::AtomID, core::Real > & atom_lddt // return by reference
) const {
	return get_stats( ref, model, res_map, true, residue_lddt, true, atom_lddt );
}

core::Real
lDDT_Calculator::get_stats(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Size > const & res_map,
	bool do_residue,
	std::map< core::Size, core::Real > & residue_stats, // return by reference
	bool do_atomistic,
	std::map< core::id::AtomID, core::Real > & atom_stats // return by reference
) const {
	lDDT_Cache cache(predicate_, ignore_oxt_);

	utility::vector1< core::Size > const & state_assignment = consider_alt_ ?
		determine_alt_states( ref, model, res_map, cache ) :
		utility::vector1< core::Size >( ref.size(), 0 ); // Just use the "zeroth" state for each entry, which is a null map.

	lDDT_Data data( do_residue, do_atomistic );
	get_stats_for_state( ref, model, res_map, state_assignment, cache, data );
	return do_division( data, residue_stats, atom_stats );
}

/// @details The approach here is to go through the full matching process, but only with the global statistics
/// Then the calling function is responsible to re-do the interation with the selected states
/// if they need the per-residue and per-atom statistics.
utility::vector1< core::Size >
lDDT_Calculator::determine_alt_states(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Size > const & res_map,
	lDDT_Cache & cache
) const {

	lDDT_Data global_data; // The data which is the same for all state combinations.

	// (Benchmarking shows that keeping track of the per-res and per-atom statistics during the state
	// finding slows things down considerably, due to the copy/summing of the residue/atom data structures.)
	// Just do the global data at first, and then come back and do the alt states later.

	// Using a flat_map instead of a regular map shaves off some time for large proteins (better memory access patterns)
	using state_data_map = boost::container::flat_map< std::pair< core::Size, core::Size >, lDDT_Data >;
	// For efficiency, we use the full pose size for the vectors, counting on most to be empty-ish.
	utility::vector1< utility::vector1< state_data_map > > res_data_map( ref.size() );

	utility::vector1< core::Size > nstate_vector( ref.size(), 1 ); // Indexed by ref residue number

	// Iteration through the res_map is the slow thing here -- move it to a flat layout to make it faster.
	// (Less jumping through memory.)
	std::vector< std::pair< core::Size, core::Size > > const flat_res_map( res_map.begin(), res_map.end() );

	for ( auto const & entry1: flat_res_map ) {
		core::Size const rres1_no = entry1.first;
		core::Size const mres1_no = entry1.second;
		core::conformation::Residue const & rres1 = ref.residue( rres1_no );

		res_data_map[ rres1_no ].resize( ref.size() );
		if ( mres1_no != 0 ) {
			nstate_vector[ rres1_no ] = cache.get_n_maps( model.residue_type( mres1_no ) );
		}

		for ( auto const & entry2: flat_res_map ) {
			core::Size const rres2_no = entry2.first;
			core::Size const mres2_no = entry2.second;

			if ( rres2_no < rres1_no ) { continue; } // Only iterate through the pairing one way (residue-internal interactions handled later)

			core::conformation::Residue const & rres2 = ref.residue( rres2_no );

			if ( rres1.polymeric_sequence_distance( rres2 ) < seqsep_ ) { continue; } // Handles self-interaction, if needed

			core::Real nbdist = rres1.nbr_atom_xyz().distance( rres2.nbr_atom_xyz() );
			if ( nbdist - rres1.nbr_radius() - rres2.nbr_radius() > R0_ ) {
				// All heavyatoms are outside the threshold - skip the atomistic enumeration
				continue;
			}

			if ( mres1_no == 0 || mres2_no == 0 ) {
				// The alternative states don't matter -- We count all the interactions as missing.
				residue_pair_stats( ref, model, rres1_no, rres2_no, mres1_no, mres2_no, global_data, cache );
				continue;
			}

			core::chemical::ResidueType const & mtype1 = model.residue_type( mres1_no );
			core::chemical::ResidueType const & mtype2 = model.residue_type( mres2_no );

			core::Size nstate1 = cache.get_n_maps( mtype1 );
			core::Size nstate2 = cache.get_n_maps( mtype2 );

			if ( nstate1 == 1 && nstate2 == 1 ) {
				// We don't have any alternate states worth mentioning
				// Just use the (default) identity map, and stick the info in the global data.
				residue_pair_stats( ref, model, rres1_no, rres2_no, mres1_no, mres2_no, global_data, cache );
				continue;
			}

			// We have multiple states. Cache the state-dependent data in the res_data_map, to be resolved later
			for ( core::Size s1(1); s1 <= nstate1; ++s1 ) {
				utility::vector1< core::Size > state1 = cache.get_mapping_state( mtype1, s1 );

				if ( mres1_no == mres2_no ) { // Same residue, so we share the state.
					lDDT_Data & data = res_data_map[ rres1_no ][ rres2_no ][ std::make_pair( s1, s1 ) ];
					residue_pair_stats( ref, model, rres1_no, rres2_no, mres1_no, mres2_no, data, cache, state1, state1 );
					continue;
				}

				for ( core::Size s2(1); s2 <= nstate2; ++s2 ) {
					utility::vector1< core::Size > state2 = cache.get_mapping_state( mtype2, s2 );

					// Will create the map if one isn't present already.
					lDDT_Data & data = res_data_map[ rres1_no ][ rres2_no ][ std::make_pair( s1, s2 ) ];
					residue_pair_stats( ref, model, rres1_no, rres2_no, mres1_no, mres2_no, data, cache, state1, state2 );

				} // state 2
			} // state 1

		} // residue 2
	} // residue 1

	///////////////////////////////////////////
	// Now we need to determine which state combination is the best one.

	// @details Doing a full combinitorial iteration of the state space is too expensive for large proteins
	// (Assuming 2 alt states for each position, N GLU/ASP/PHE/TYR/ARG residues would result in 2^N possibilities, which rapidly becomes too much.)
	// Instead, we assume we can simply hill-climb it: Each state swap would result in

	utility::vector1< core::Size > state_assignment( ref.size(), 1 ); // Indexed by ref residue number (first time through will replace all.
	debug_assert( nstate_vector.size() == state_assignment.size() );

	// Get the info for the current best state.
	lDDT_Data best_state;
	get_stats_for_state( ref, model, res_map, state_assignment, cache, best_state );
	core::Real best_global = best_state.global_lddt();

	bool found_improvement = true;

	while ( found_improvement ) {
		found_improvement = false;

		for ( core::Size ii(1); ii <= state_assignment.size(); ++ii ) {
			for ( core::Size ss(1); ss <= nstate_vector[ii]; ++ss ) {
				core::Size const old_ss = state_assignment[ii];
				if ( ss == old_ss ) { continue; }

				lDDT_Data current_data( best_state ); // Make a copy of the intermediate states.

				// Now we need to a) remove the data for the old state and b) add the data for the new state.

				core::Size const rres1_no = ii; // If we have states to swap, we're in the map.

				for ( auto const & entry2: flat_res_map ) {
					core::Size const rres2_no = entry2.first;

					auto const & res1_data_map = res_data_map[ std::min(rres1_no,rres2_no) ];
					if ( res1_data_map.empty() ) {
						continue; // No alt states for this
					}
					auto const & res12_data_map = res1_data_map[ std::max(rres1_no,rres2_no) ];
					if ( res12_data_map.empty() ) {
						continue; // No alt states for this.
					}

					core::Size other_state = state_assignment[ rres2_no ];

					auto old_state_key = (rres1_no < rres2_no) ? std::make_pair( old_ss, other_state ) : std::make_pair( other_state, old_ss );
					auto new_state_key = (rres1_no < rres2_no) ? std::make_pair( ss, other_state ) : std::make_pair( other_state, ss );

					if ( rres1_no == rres2_no ) {
						old_state_key = std::make_pair( old_ss, old_ss );
						new_state_key = std::make_pair( ss, ss );
					}

					auto const & old_iter = res12_data_map.find( old_state_key );
					if ( old_iter != res12_data_map.end() ) {
						current_data -= old_iter->second;
					}

					auto const & new_iter = res12_data_map.find( new_state_key );
					if ( new_iter != res12_data_map.end() ) {
						current_data += new_iter->second;
					}
				}

				// The data should be updated, so check it for improvement
				core::Real const & new_global = current_data.global_lddt();
				if ( new_global > best_global ) {
					found_improvement = true; // Do another cycle of checking.
					best_global = new_global;
					state_assignment[ ii ] = ss;
					best_state = std::move( current_data ); // Don't need the old one.
				}
			}
		}
	}

	return state_assignment;
}

core::Real
lDDT_Calculator::get_stats_for_state(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	std::map< core::Size, core::Size > const & res_map,
	utility::vector1< core::Size > const & state_assignment,
	lDDT_Cache & cache,
	lDDT_Data & data
) const {

	// Iterate through all the atoms in the reference pose, looking for pairs within the threshold.
	for ( auto const & entry1: res_map ) {
		core::Size const rres1_no = entry1.first;
		core::Size const mres1_no = entry1.second;

		core::conformation::Residue const & rres1 = ref.residue( rres1_no );

		utility::vector1< core::Size > const & mres1_map =
			( mres1_no != 0 )
			? cache.get_mapping_state( model.residue_type( mres1_no), state_assignment[ mres1_no ] )
			: lDDT_Cache::identity_map_;

		for ( auto const & entry2: res_map ) {
			core::Size const rres2_no = entry2.first;
			core::Size const mres2_no = entry2.second;

			if ( rres2_no < rres1_no ) { continue; } // Only iterate through the pairing one way (residue-internal interactions handled later)

			core::conformation::Residue const & rres2 = ref.residue( rres2_no );

			if ( rres1.polymeric_sequence_distance( rres2 ) < seqsep_ ) { continue; } // Handles self-interaction, if needed

			core::Real nbdist = rres1.nbr_atom_xyz().distance( rres2.nbr_atom_xyz() );
			if ( nbdist - rres1.nbr_radius() - rres2.nbr_radius() > R0_ ) {
				// All heavyatoms are outside the threshold - skip the atomistic enumeration
				continue;
			}

			utility::vector1< core::Size > const & mres2_map =
				( mres2_no != 0 )
				? cache.get_mapping_state( model.residue_type( mres2_no), state_assignment[ mres2_no ] )
				: lDDT_Cache::identity_map_;

			// Assemble the statistics for this residue pair into data.
			residue_pair_stats( ref, model, rres1_no, rres2_no, mres1_no, mres2_no, data, cache, mres1_map, mres2_map );

		} // for residue 2
	} // for residue 1

	if ( data.global_ndist ) {
		return core::Real( data.global_nmatch ) / data.global_ndist;
	} else {
		return 0.0; // Avoid division by zero
	}
}


core::Real
lDDT_Calculator::do_division(
	lDDT_Data & data,
	std::map< core::Size, core::Real > & residue_stats,
	std::map< core::id::AtomID, core::Real > & atom_stats
) const {
	// We should only get entries in the data maps if we have non-zero counts.

	residue_stats.clear();
	for ( auto const & entry: data.res_nmatch ) {
		TR.Trace << "RESIDUE " << entry.first << " Conserved: " << entry.second << " out of " << data.res_ndist[ entry.first ] << std::endl;
		residue_stats[ entry.first ] = core::Real( entry.second )/data.res_ndist[ entry.first ];
	}

	atom_stats.clear();
	for ( auto const & entry: data.atom_nmatch ) {
		atom_stats[ entry.first ] = core::Real( entry.second )/data.atom_ndist[ entry.first ];
	}

	return data.global_lddt();
}

void
lDDT_Calculator::residue_pair_stats(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	core::Size rres1_no,
	core::Size rres2_no,
	core::Size mres1_no,
	core::Size mres2_no,
	lDDT_Data & data,
	lDDT_Cache & cache,
	utility::vector1< core::Size > const & mres1_map,
	utility::vector1< core::Size > const & mres2_map
) const {

	core::Size const nthresholds = thresholds_.size();
	core::conformation::Residue const & rres1 = ref.residue( rres1_no );
	core::conformation::Residue const & rres2 = ref.residue( rres2_no );

	for ( core::Size ratm1: cache.atoms_matching_predicate( ref, rres1_no ) ) {
		for ( core::Size ratm2: cache.atoms_matching_predicate( ref, rres2_no ) ) {

			if ( rres1_no == rres2_no && ratm2 <= ratm1 ) { continue; } // Handle within-residue stats, if necessary

			core::Real const refdist = rres1.xyz( ratm1 ).distance( rres2.xyz( ratm2 ) );
			if ( refdist > R0_ ) { continue; }

			// Calculate the delta value (-1 means that we're an invalid distance on the model)
			core::Real delta = -1;
			if ( mres1_no != 0 && mres2_no != 0 ) {
				core::conformation::Residue const & mres1 = model.residue( mres1_no );
				core::conformation::Residue const & mres2 = model.residue( mres2_no );

				core::Size matm1 = get_matching_atom(rres1, mres1, ratm1, mres1_map);
				core::Size matm2 = get_matching_atom(rres2, mres2, ratm2, mres2_map);

				if ( matm1 != 0 && matm2 != 0 ) {
					core::Real const model_dist = mres1.xyz( matm1 ).distance( mres2.xyz( matm2 ) );
					delta = std::abs( refdist - model_dist );
				}
			}

			// Handle incrementing counts;
			data.global_ndist += nthresholds;
			if ( data.do_residue ) {
				data.res_ndist[ rres1_no ] += nthresholds; // defaults zero
				data.res_ndist[ rres2_no ] += nthresholds;
			}
			if ( data.do_atomistic ) {
				data.atom_ndist[ id::AtomID( ratm1, rres1_no) ] += nthresholds;
				data.atom_ndist[ id::AtomID( ratm2, rres2_no) ] += nthresholds;
			}

			if ( delta < 0 ) { continue; }

			for ( core::Real const threshold: thresholds_ ) {
				if ( delta > threshold ) {
					continue;
				}

				++data.global_nmatch;

				if ( data.do_residue ) {
					++data.res_nmatch[ rres1_no ];
					++data.res_nmatch[ rres2_no ];
				}

				if ( data.do_atomistic ) {
					++data.atom_nmatch[ id::AtomID( ratm1, rres1_no) ];
					++data.atom_nmatch[ id::AtomID( ratm2, rres2_no) ];
				}
			} // for thresholds
		} // for atom2
	} // for atom 1
}


std::map< core::Size, core::Size >
lDDT_Calculator::make_identity_map( core::pose::Pose const & pose ) const {
	std::map< core::Size, core::Size > mapping;
	for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
		if ( pose.residue_type(ii).is_virtual_residue() ) { continue; } // Ignore virtual residues by default.
		mapping[ ii ] = ii;
	}
	return mapping;
}

core::Size
lDDT_Calculator::get_matching_atom(
	core::conformation::Residue const & rres,
	core::conformation::Residue const & mres,
	core::Size ratm,
	utility::vector1< core::Size > const & mres_map
) const {
	if ( ratm == 0 ) { return 0; }

	core::Size matm = ratm; // Assume easy case of same residue type
	if ( &rres.type() != &mres.type() ) { // Address-based comparison for speed.
		// Do harder name-based approach
		std::string const & name = rres.atom_name( ratm );
		if ( ! mres.has( name ) ) { return 0; }
		matm = mres.atom_index( name );
	}

	if ( matm == 0 || matm > mres_map.size() ) {
		return matm;
	} else {
		return mres_map[ matm ];
	}
}

/////////////////  lDDT_Cache ///////////////////////////

// Declaration of static data member.
utility::vector1< core::Size > const lDDT_Calculator::lDDT_Cache::identity_map_;


utility::vector1< core::Size > const &
lDDT_Calculator::lDDT_Cache::atoms_matching_predicate(
	core::pose::Pose const & pose,
	core::Size resno
) {
	core::conformation::Residue const * ptr = & pose.residue(resno);

	if ( predicate_atoms_.count(ptr) == 0 ) {
		// Fill cached values
		utility::vector1< core::Size > & atoms = predicate_atoms_[ ptr ]; // Reference to in-place
		for ( core::Size atm(1); atm <= ptr->natoms(); ++atm ) {
			if ( ptr->is_virtual( atm ) ) { continue; } // Skip virtual atoms
			if ( ignore_oxt_ && ptr->has( "OXT" ) && ptr->atom_index( "OXT" ) == atm ) { continue; }
			// Technically should take two different poses, but all the actual ones don't look at the second.
			if ( (*predicate_)(pose, pose, resno, atm ) ) {
				atoms.push_back( atm );
			}
		}
	}

	return predicate_atoms_[ ptr ];
}

core::Size
lDDT_Calculator::lDDT_Cache::get_n_maps(
	core::chemical::ResidueType const & type
) {
	auto const & states = residue_mapping_states(type);
	if ( states.empty() ) {
		return 1;
	} else {
		return states.size();
	}
}

utility::vector1< utility::vector1< core::Size > > const &
lDDT_Calculator::lDDT_Cache::residue_mapping_states(
	core::chemical::ResidueType const & type
) {
	if ( mapping_states_.count( &type ) == 0 ) {
		make_mapping_states(type);
	}
	return mapping_states_[ &type ];
}

utility::vector1< core::Size > const &
lDDT_Calculator::lDDT_Cache::get_mapping_state(
	core::chemical::ResidueType const & type,
	core::Size state_no
) {
	if ( state_no == 0 ) {
		return identity_map_;
	}
	utility::vector1< utility::vector1< core::Size > > const & all_states = residue_mapping_states( type );
	if ( state_no > all_states.size() ) {
		return identity_map_;
	}
	return all_states[ state_no ];
}

bool
lDDT_Calculator::lDDT_Cache::swap_indexes(
	utility::vector1< core::Size > & vec,
	core::chemical::ResidueType const & type,
	char const * const atm1, // String Literal
	char const * const atm2
) {
	if ( type.has( atm1 ) && type.has( atm2 ) ) {
		vec[ type.atom_index(atm1) ] = type.atom_index(atm2);
		vec[ type.atom_index(atm2) ] = type.atom_index(atm1);
		return true;
	} else {
		return false;
	}
}

void
lDDT_Calculator::lDDT_Cache::add_states_from_automorphisms(
	utility::vector1< utility::vector1< core::Size > > & states,
	core::chemical::ResidueType const & type
) {
	states.clear(); // We're going to replace the implict version

	core::chemical::AutomorphismIterator iter( type, /*includeH*/ false );

	while ( true ) {
		states.emplace_back( iter.next() );
		if ( states.back().size() == 0 ) {
			// Empty vector indicates end of iteration.
			states.pop_back();
			break;
		}
	}

	TR << "ResidueType " << type.name() << " has " << states.size() << " automorphic states." << std::endl;
}

void
lDDT_Calculator::lDDT_Cache::make_mapping_states(
	core::chemical::ResidueType const & type
) {
	using namespace core::chemical;

	utility::vector1< utility::vector1< core::Size > > & states = mapping_states_[ &type ];

	// First deal with the unity state
	states.emplace_back( utility::vector1< core::Size >{} );
	utility::vector1< core::Size > & state1 = states[1];
	for ( core::Size ii(1); ii <= type.nheavyatoms(); ++ii ) {
		state1.push_back(ii);
	}

	// Deal with the alternate states

	// Somewhat hacky approach -- attempt to hard code the standard types.
	switch ( type.aa() ) {
	// Most of the standard AAs don't have any equivalents other than the atom itself.
	case aa_ala:
	case aa_cys:
	case aa_gly:
	case aa_his:
	case aa_ile:
	case aa_lys:
	case aa_met:
	case aa_asn:
	case aa_pro:
	case aa_gln:
	case aa_ser:
	case aa_thr:
	case aa_trp:
	case aa_dal : // D
	case aa_dcs :
		// No D-Gly
	case aa_dhi:
	case aa_dil:
	case aa_dly:
	case aa_dme:
	case aa_dan:
	case aa_dpr:
	case aa_dgn:
	case aa_dse:
	case aa_dth:
	case aa_dtr:
	case aa_b3a : // Beta
	case aa_b3c:
	case aa_b3g:
	case aa_b3h:
	case aa_b3i:
	case aa_b3k:
	case aa_b3m:
	case aa_b3n:
	case aa_b3p:
	case aa_b3q:
	case aa_b3s:
	case aa_b3t:
	case aa_b3w:
	case ou3_ala : // Oligourea
	case ou3_cys:
	case ou3_gly:
	case ou3_his:
	case ou3_ile:
	case ou3_lys:
	case ou3_met:
	case ou3_asn:
	case ou3_pro:
	case ou3_gln:
	case ou3_ser:
	case ou3_thr:
	case ou3_trp :
		break;
		// The original paper indicates that val & leu have chemically equivalent atoms,
		// but that's not really the case - while they're not chiral, they are pro-chiral.
		// The Swiss-Model online calculator seems to (correctly) ignore VL alternate states.
	case aa_val:
	case aa_leu:
	case aa_dva : // Others
	case aa_dle:
	case aa_b3v:
	case aa_b3l:
	case ou3_val:
	case ou3_leu :
		break;
		// Phe & Tyr - ring flip
	case aa_phe:
	case aa_tyr:
	case aa_dph : // Others
	case aa_dty:
	case aa_b3f:
	case aa_b3y:
	case ou3_phe:
	case ou3_tyr :
		states.emplace_back( state1 );
		if ( !( swap_indexes( states.back(), type, "CD1", "CD2" ) &&
				swap_indexes( states.back(), type, "CE1", "CE2" )
				) ) {
			states.pop_back(); // clean up error
		}
		break;
		// Asp - carboxylate flip
	case aa_asp:
	case aa_das:
	case aa_b3d:
	case ou3_asp :
		states.emplace_back( state1 );
		if ( !swap_indexes( states.back(), type, "OD1", "OD2" ) ) {
			states.pop_back(); // clean up error
		}
		break;
		// Glu - carboxylate flip
	case aa_glu:
	case aa_dgu:
	case aa_b3e:
	case ou3_glu :
		states.emplace_back( state1 );
		if ( !swap_indexes( states.back(), type, "OE1", "OE2" ) ) {
			states.pop_back(); // clean up error
		}
		break;
		// Arg - guanidinium flip
	case aa_arg:
	case aa_dar:
	case aa_b3r:
	case ou3_arg :
		states.emplace_back( state1 );
		if ( !swap_indexes( states.back(), type, "NH1", "NH22" ) ) {
			states.pop_back(); // clean up error
		}
		break;
		// Nucleic acids -- The oxygens on the phosphate are chemically equivalent
	case na_ade:
	case na_cyt:
	case na_gua:
	case na_thy:
	case na_rad:
	case na_rcy:
	case na_rgu:
	case na_ura:
	case na_lra:
	case na_lrc:
	case na_lrg:
	case na_lur :
		states.emplace_back( state1 );
		if ( !swap_indexes( states.back(), type, "OP1", "OP2" ) ) {
			states.pop_back(); // clean up error
		}
		break;
		// Ligands -- handle them specially.
	case aa_unk :
		if ( type.is_ligand() ) {
			add_states_from_automorphisms( states, type );
		}
		break;
		// Not otherwise specified? Just skip the equivalences
	default :
		break;
	};

	// If we're not ignoring OXT, account for the OXT/O swap
	if ( (!ignore_oxt_) && type.has("OXT") && type.has("O") ) {
		// Make a copy of the states vector to avoid modifying what we're iterating over.
		utility::vector1< utility::vector1< core::Size > > states_copy( states );
		for ( auto const & state: states_copy ) {
			states.emplace_back( state );
			if ( !swap_indexes( states.back(), type, "OXT", "O" ) ) {
				states.pop_back();
			}
		}
	}

}

/////////////////  lDDT_Data ///////////////////////////

core::Real
lDDT_Calculator::lDDT_Data::global_lddt() const {
	if ( global_ndist ) {
		return core::Real( global_nmatch ) / global_ndist;
	} else {
		return 0.0; // Avoid division by zero
	}
}

template< typename T >
void
lDDT_Calculator::lDDT_Data::add_inplace( std::map< T, core::Size> & lhs, std::map< T, core::Size > const & rhs ) {
	for ( auto const & entry: rhs ) {
		lhs[ entry.first ] += entry.second; // Will zero initialize, if necessary.
	}
}

template< typename T >
void
lDDT_Calculator::lDDT_Data::subtract_inplace( std::map< T, core::Size> & lhs, std::map< T, core::Size > const & rhs ) {
	for ( auto const & entry: rhs ) {
		lhs[ entry.first ] -= entry.second;
	}
}

void
lDDT_Calculator::lDDT_Data::operator+=( lDDT_Data const & rhs ) {
	global_ndist += rhs.global_ndist;
	global_nmatch += rhs.global_nmatch;

	add_inplace( res_ndist, rhs.res_ndist );
	add_inplace( res_nmatch, rhs.res_nmatch );

	add_inplace( atom_ndist, rhs.atom_ndist );
	add_inplace( atom_nmatch, rhs.atom_nmatch );
}

void
lDDT_Calculator::lDDT_Data::operator-=( lDDT_Data const & rhs ) {
	global_ndist -= rhs.global_ndist;
	global_nmatch -= rhs.global_nmatch;

	subtract_inplace( res_ndist, rhs.res_ndist );
	subtract_inplace( res_nmatch, rhs.res_nmatch );

	subtract_inplace( atom_ndist, rhs.atom_ndist );
	subtract_inplace( atom_nmatch, rhs.atom_nmatch );
}

} // namespace scoring
} // namespace core
