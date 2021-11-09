// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/util.cc
/// @brief  Nonmember functions for evaluating some or all energy methods on residues or residue pairs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/util.hh>

// Package headers
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>


// Numeric headers

static basic::Tracer TR( "core.scoring.utils" );

namespace core {
namespace scoring {

void
eval_scsc_sr2b_energies(
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	Vector const & r1sc_centroid,
	Vector const & r2sc_centroid,
	Real const & r1sc_radius,
	Real const & r2sc_radius,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) {
	Real const scsc_d2 = r1sc_centroid.distance_squared( r2sc_centroid );
	Real const scsc_radsum = r1sc_radius + r2sc_radius;
	for ( auto
			iter = sfxn.ci_2b_begin(), iter_end = sfxn.ci_2b_end();
			iter != iter_end; ++iter ) {
		Real cutoff = scsc_radsum + (*iter)->atomic_interaction_cutoff();
		if ( ! (*iter)->divides_backbone_and_sidechain_energetics() || scsc_d2 < cutoff * cutoff ) {
			(*iter)->sidechain_sidechain_energy( r1, r2, pose, sfxn, emap );
		}
	}
	for ( auto
			iter = sfxn.cd_2b_begin(), iter_end = sfxn.cd_2b_end();
			iter != iter_end; ++iter ) {
		Real cutoff = scsc_radsum + (*iter)->atomic_interaction_cutoff();
		if ( ! (*iter)->divides_backbone_and_sidechain_energetics() || scsc_d2 < cutoff * cutoff ) {
			(*iter)->sidechain_sidechain_energy( r1, r2, pose, sfxn, emap );
		}
	}
}

void
eval_bbsc_sr2b_energies(
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	Vector const & r1bb_centroid,
	Vector const & r2sc_centroid,
	Real const & r1bb_radius,
	Real const & r2sc_radius,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) {
	Real const bbsc_d2 = r1bb_centroid.distance_squared( r2sc_centroid );
	Real const bbsc_radsum = r1bb_radius + r2sc_radius;
	for ( auto
			iter = sfxn.ci_2b_begin(), iter_end = sfxn.ci_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbsc_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbsc_d2 < cutoff * cutoff ) {
				(*iter)->backbone_sidechain_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
	for ( auto
			iter = sfxn.cd_2b_begin(), iter_end = sfxn.cd_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbsc_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbsc_d2 < cutoff * cutoff ) {
				(*iter)->backbone_sidechain_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
}

void
eval_bbbb_sr2b_energies(
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	Vector const & r1bb_centroid,
	Vector const & r2bb_centroid,
	Real const & r1bb_radius,
	Real const & r2bb_radius,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) {
	Real const bbbb_d2 = r1bb_centroid.distance_squared( r2bb_centroid );
	Real const bbbb_radsum = r1bb_radius + r2bb_radius;
	for ( auto
			iter = sfxn.ci_2b_begin(), iter_end = sfxn.ci_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbbb_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbbb_d2 < cutoff * cutoff ) {
				(*iter)->backbone_backbone_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
	for ( auto
			iter = sfxn.cd_2b_begin(), iter_end = sfxn.cd_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbbb_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbbb_d2 < cutoff * cutoff ) {
				(*iter)->backbone_backbone_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
}

/// @details returns the origin if there are no backbone atoms
Vector
compute_bb_centroid(
	conformation::Residue const & res
) {
	Vector bb_centroid( 0.0 );
	Size count_n_bb( 0 );
	for ( Size ii = 1; ii <= res.type().first_sidechain_atom() - 1; ++ii ) {
		++count_n_bb;
		bb_centroid += res.xyz( ii );
	}
	if ( count_n_bb != 0 ) bb_centroid /= count_n_bb;
	return bb_centroid;
}

Real
compute_bb_radius(
	conformation::Residue const & res,
	Vector const & bb_centroid
) {
	Real bb_radius = 0;
	for ( Size ii = 1; ii <= res.type().first_sidechain_atom() - 1; ++ii ) {
		Real d2 = res.xyz( ii ).distance_squared( bb_centroid );
		if ( bb_radius < d2 ) bb_radius = d2;
	}
	return bb_radius;
}

Vector
compute_sc_centroid(
	conformation::Residue const & res
) {
	Vector centroid( 0.0 );
	Size count( 0 );
	for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
		count += 1;
		centroid += res.xyz( ii );
	}
	if ( count == 0 ) {
		return compute_bb_centroid( res );
	} else {
		centroid /= count;
		return centroid;
	}
}

Real
compute_sc_radius(
	conformation::Residue const & res,
	Vector const & centroid
) {
	Real max_d2 = 0;
	for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
		Real d2 = res.xyz(ii).distance_squared( centroid );
		if ( d2 > max_d2 ) max_d2 = d2;
	}
	return std::sqrt( max_d2 );
}

bool
check_score_function_sanity(
	utility::options::OptionCollection const & options,
	std::string const & scorefxn_key,
	bool throw_exception
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( utility::startswith( scorefxn_key, "score12" ) || utility::startswith( scorefxn_key, "pre_talaris" ) ) {
		if ( ! options[ mistakes::restore_pre_talaris_2013_behavior ] ||
				options[ corrections::restore_talaris_behavior ] ) {
			TR.Error
				<< "**********************************************" << std::endl
				<< "To use the '" << scorefxn_key << "' score function" << std::endl
				<< "please use the " << std::endl
				<< std::endl
				<< "        -restore_pre_talaris_2013_behavior" << std::endl
				<< std::endl
				<< "flag on the command line." << std::endl
				<< std::endl
				<< "(And do not use the -restore_talaris_behavior flag.)" << std::endl
				<< "**********************************************" << std::endl;
			if ( throw_exception ) {
				throw CREATE_EXCEPTION(utility::excn::BadInput, "Missing -restore_pre_talaris_2013_behavior flag.");
			} else {
				return false;
			}
		}
	} else if ( utility::startswith( scorefxn_key, "talaris20" ) ) {
		if ( options[ mistakes::restore_pre_talaris_2013_behavior ] ||
				! options[ corrections::restore_talaris_behavior ] ) {
			TR.Error
				<< "**********************************************" << std::endl
				<< "To use the '" << scorefxn_key << "' score function" << std::endl
				<< "please use the " << std::endl
				<< std::endl
				<< "        -restore_talaris_behavior" << std::endl
				<< std::endl
				<< "flag on the command line." << std::endl
				<< std::endl
				<< "(And do not use the -restore_pre_talaris_2013_behavior flag.)" << std::endl
				<< "**********************************************" << std::endl;
			if ( throw_exception ) {
				throw CREATE_EXCEPTION(utility::excn::BadInput, "Missing -restore_talaris_behavior flag.");
			} else {
				return false;
			}
		}
	} else if ( utility::startswith( scorefxn_key, "ref20" ) || utility::startswith( scorefxn_key, "beta_nov15" ) ) {
		if ( options[ mistakes::restore_pre_talaris_2013_behavior ] ||
				options[ corrections::restore_talaris_behavior ] ) {
			TR.Error
				<< "**********************************************" << std::endl
				<< "To use the '" << scorefxn_key << "' score function" << std::endl
				<< "please do not use either of the" << std::endl
				<< std::endl
				<< "        -restore_talaris_behavior" << std::endl
				<< "        -restore_pre_talaris_2013_behavior" << std::endl
				<< std::endl
				<< "flags on the command line." << std::endl
				<< "**********************************************" << std::endl;
			if ( throw_exception ) {
				throw CREATE_EXCEPTION(utility::excn::BadInput, "Using -restore_talaris_behavior or -restore_pre_talaris_2013_behavior with REF-era scorefunctions.");
			} else {
				return false;
			}
		}
	}
	// Likely sane score function
	return true;
}


std::map< id::AtomID, utility::vector1< Real > >
get_single_atom_energies(
	core::pose::Pose & pose,
	ScoreFunction const & scorefxn,
	ScoreTypes const & types,
	select::residue_selector::ResidueSubsetOP subset
) {
	if ( subset != nullptr && subset->size() != pose.size() ) {
		utility_exit_with_message("Ill-formed residue subset: " + std::to_string( subset->size() ) + " entries for a pose of size " + std::to_string( pose.size() ) );
	}

	std::map< id::AtomID, utility::vector1< Real > > energy_map;

	scorefxn( pose ); // Rescore to make sure we have proper scoring setup.

	core::scoring::EnergyMap const & weights( scorefxn.weights() );
	utility::vector1< methods::EnergyMethodCOP > onebody_terms;
	for ( auto iter = scorefxn.all_energies_begin(); iter != scorefxn.all_energies_end(); ++iter ) {
		if ( (*iter)->has_atomistic_energies() ) {
			onebody_terms.push_back( *iter );
		}
	}

	if ( onebody_terms.empty() ) { return energy_map; } // Nothing to do.

	for ( core::Size ii(1); ii <= pose.size(); ++ii ) { // Residues
		if ( subset && ! (*subset)[ii] ) { continue; } // Ignore if not part of selection

		core::conformation::Residue const & res( pose.residue(ii) );
		for ( core::Size aa(1); aa <= res.natoms(); ++aa ) {
			EnergyMap onebody; // It accumulates
			for ( auto const & term: onebody_terms ) {
				term->atomistic_energy( aa, res, pose, scorefxn, onebody );
			}
			onebody *= weights;
			if ( onebody.sum() == 0.0 ) { continue; }

			utility::vector1< Real > energies;
			for ( ScoreType st: types ) {
				energies.push_back( onebody[ st ] );
			}

			energy_map[ id::AtomID( aa, ii ) ] = energies;
		}
	}

	return energy_map;
}

std::map< std::pair< id::AtomID, id::AtomID >,  utility::vector1< Real > >
get_pairwise_atom_energies(
	core::pose::Pose & pose,
	ScoreFunction const & scorefxn,
	ScoreTypes const & types,
	select::residue_selector::ResidueSubsetOP subset1,
	select::residue_selector::ResidueSubsetOP subset2
) {
	if ( subset1 != nullptr && subset1->size() != pose.size() ) {
		utility_exit_with_message("Ill-formed residue subset: " + std::to_string( subset1->size() ) + " entries for a pose of size " + std::to_string( pose.size() ) );
	}
	if ( subset2 != nullptr && subset2->size() != pose.size() ) {
		utility_exit_with_message("Ill-formed secondary residue subset: " + std::to_string( subset2->size() ) + " entries for a pose of size " + std::to_string( pose.size() ) );
	}

	std::map< std::pair< id::AtomID, id::AtomID >, utility::vector1< Real > > energy_map;

	scorefxn( pose ); // Rescore to make sure we have proper scoring setup.

	core::scoring::EnergyMap const & weights( scorefxn.weights() );
	utility::vector1< methods::EnergyMethodCOP > pairwise_terms;
	for ( auto iter = scorefxn.all_energies_begin(); iter != scorefxn.all_energies_end(); ++iter ) {
		if ( (*iter)->has_atomistic_pairwise_energies() ) {
			pairwise_terms.push_back( *iter );
		}
	}

	if ( pairwise_terms.empty() ) { return energy_map; } // Nothing to do.

	for ( core::Size ii(1); ii <= pose.size(); ++ii ) { // Residues
		// At this point, only ignore if we've defined both subsets and it's in neither
		if ( subset1 && subset2 && !(*subset1)[ii] && !(*subset2)[ii] ) { continue; }

		core::conformation::Residue const & res1( pose.residue(ii) );
		for ( core::Size jj(ii); jj <= pose.size(); ++jj ) { // Residues (including the self interaction
			// If a subset is given, at least one partner must be in it.
			if ( subset1 && !(*subset1)[ii] && !(*subset1)[jj] ) { continue; }
			if ( subset2 && !(*subset2)[ii] && !(*subset2)[jj] ) { continue; }

			core::conformation::Residue const & res2( pose.residue(jj) );
			for ( core::Size aa(1); aa <= res1.natoms(); ++aa ) {
				for ( core::Size bb(1); bb <= res2.natoms(); ++bb ) {
					if ( ii == jj && bb <= aa ) { continue; } // Only consider the same-residue pairing once.

					EnergyMap twobody; // It accumulates
					for ( auto const & term: pairwise_terms ) {
						term->atomistic_pair_energy( aa, res1, bb, res2, pose, scorefxn, twobody );
					}

					twobody *= weights;
					if ( twobody.sum() == 0.0 ) { continue; } // Don't include pairs with zero energy.

					utility::vector1< Real > energies;
					for ( ScoreType st: types ) {
						energies.push_back( twobody[ st ] );
					}

					energy_map[ std::make_pair( id::AtomID( aa, ii ), id::AtomID( bb, jj ) ) ] = energies;
				} // bb
			} // aa
		} // jj
	} // ii

	return energy_map;
}

std::map< id::AtomID, utility::vector1< Real > >
get_atomistic_energies(
	core::pose::Pose & pose,
	ScoreFunction const & scorefxn,
	ScoreTypes const & types,
	select::residue_selector::ResidueSubsetOP subset
) {
	std::map< id::AtomID,  utility::vector1< Real > > energy_map = get_single_atom_energies(pose, scorefxn, types, subset );

	for ( auto const & entry: get_pairwise_atom_energies(pose, scorefxn, types, subset ) ) {
		id::AtomID const & id1 = entry.first.first;
		id::AtomID const & id2 = entry.first.second;
		utility::vector1< Real > energies = entry.second;
		// Split the energies between the two atoms in the pair
		for ( Size ii(1); ii <= energies.size(); ++ii ) {
			energies[ii] /= 2.0;
		}

		if ( energy_map.count( id1 ) == 0 ) {
			energy_map[ id1 ] = energies;
		} else {
			for ( Size ii(1); ii <= energies.size(); ++ii ) {
				energy_map[ id1 ][ii] += energies[ii];
			}
		}

		if ( energy_map.count( id2 ) == 0 ) {
			energy_map[ id2 ] = energies;
		} else {
			for ( Size ii(1); ii <= energies.size(); ++ii ) {
				energy_map[ id2 ][ii] += energies[ii];
			}
		}
	}

	return energy_map;
}


}
}
