// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/carbohydrates/util.cc
/// @brief   Utility function definitions for scoring carbohydrate-containing poses.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit Headers
#include <core/scoring/carbohydrates/util.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/exit.hh>


namespace core {
namespace scoring {
namespace carbohydrates {
using namespace core::conformation::carbohydrates;

// Get the CHI Energy Function linkage type for phi for a particular residue.
/// @return  CHIEnergyFunctionLinkageType, which corresponds to one of two particular CHI Energy Functions specific to
/// phi, or LINKAGE_NA if the CHI Energy Function is not applicable to this torsion.
CHIEnergyFunctionLinkageType
get_CHI_energy_function_linkage_type_for_phi_for_residue_in_pose( pose::Pose const & pose, core::uint rsd_num )
{
	conformation::Residue const & rsd( pose.residue( rsd_num ) );

	if ( rsd.is_carbohydrate() && rsd.seqpos() != 1 ) {
		chemical::carbohydrates::CarbohydrateInfoCOP info( rsd.carbohydrate_info() );
		if ( info->is_pyranose() ) {
			// TODO: Wood's lab assumed that the rings would always be 4C1 chairs. If an alpha sugar is flipped, it
			// should probably be treated as a beta.  I should probably abandon getting the anomeric form and explicitly
			// determine axial/equatorial here too.
			if ( info->is_alpha_sugar() ) {
				return ALPHA_LINKS;
			} else if ( info->is_beta_sugar() ) {
				return BETA_LINKS;
			}
		}
	}
	return LINKAGE_NA;
}

// Get the CHI Energy Function linkage type for phi for a particular residue.
/// @return  CHIEnergyFunctionLinkageType, which corresponds to one of four particular CHI Energy Functions specific to
/// psi, or LINKAGE_NA if the CHI Energy Function is not applicable to this torsion.
CHIEnergyFunctionLinkageType
get_CHI_energy_function_linkage_type_for_psi_for_residue_in_pose( pose::Pose const & pose, core::uint rsd_num )
{
	using namespace chemical::rings;
	using namespace conformation;

	Residue const & rsd( pose.residue( rsd_num ) );

	if ( rsd.is_carbohydrate() && rsd.seqpos() != 1 ) {
		// For psi, we need to get information from the previous residue.

		core::Size prev_rsd_num;
		if ( pose.glycan_tree_set() ) {
			prev_rsd_num = pose.glycan_tree_set()->get_parent( rsd_num ) ;
		} else {
			prev_rsd_num = find_seqpos_of_saccharides_parent_residue( rsd );
		}

		// If this is not a saccharide residue, do nothing.

		Residue const & prev_rsd( pose.residue( prev_rsd_num ));
		if ( ! prev_rsd.is_carbohydrate() ) { return LINKAGE_NA; }

		if ( prev_rsd.type().is_cyclic() ) {
			// What is our connecting atom?
			uint const connect_atom( prev_rsd.connect_atom( rsd ) );

			// Now, get the ring atoms.
			// We can assume that the ring we care about is the 1st ring, since this is a sugar.
			utility::vector1< uint > const ring_atoms( prev_rsd.type().ring_atoms( 1 ) );

			// Next, we must figure out which position on the ring has the glycosidic bond.
			uint const position( position_of_atom_on_ring( prev_rsd, connect_atom, ring_atoms ) );
			if ( position != 0 ) {  // We found a position, so this is not an exocyclic linkage.
				// If this is not a pyranose, do nothing, because we do not have QM data for this residue.
				if ( ! prev_rsd.carbohydrate_info()->is_pyranose() ) { return LINKAGE_NA; }

				// TODO: This needs to shift odd and even for n-ketoses where n is even.
				if ( prev_rsd.carbohydrate_info()->is_ketose() ) { return LINKAGE_NA; }  // TEMP

				// Finally, check if it's axial or equatorial and call the appropriate function.
				switch ( is_atom_axial_or_equatorial_to_ring( prev_rsd, connect_atom, ring_atoms ) ) {
				case AXIAL :
					if ( position % 2 == 0 ) {  // even
						return _2AX_3EQ_4AX_LINKS;
					} else /* odd */ {
						return _2EQ_3AX_4EQ_LINKS;
					}
					break;
				case EQUATORIAL :
					if ( position % 2 == 0 ) {  // even
						return _2EQ_3AX_4EQ_LINKS;
					} else /* odd */ {
						return _2AX_3EQ_4AX_LINKS;
					}
					break;
				case NEITHER :
					return LINKAGE_NA;
				}
			}
		}  // prev_rsd.type().is_cyclic()

		// If we've reached this point, we are either connected to an acyclic sugar or we are at an exocyclic linkage.
		// In such cases, psi's preferences depend on whether we are alpha or beta.
		chemical::carbohydrates::CarbohydrateInfoCOP info( rsd.carbohydrate_info() );
		if ( info->is_pyranose() ) {
			if ( info->is_alpha_sugar() ) {
				return ALPHA6_LINKS;
			} else if ( info->is_beta_sugar() ) {
				return BETA6_LINKS;
			}
		}
	}
	return LINKAGE_NA;
}

// Get the omega preference for a particular residue.
OmegaPreferenceType
get_omega_preference_for_residue_in_pose( pose::Pose const & pose, core::uint rsd_num )
{
	using namespace chemical::rings;
	using namespace conformation;

	Residue const & rsd( pose.residue( rsd_num ) );

	//JAB - shouldn't we check parent residue here and make sure its not zero instead of 1?
	if ( rsd.is_carbohydrate() && rsd.seqpos() != 1 ) {
		// For omega, we need to get information from the previous residue.

		core::Size prev_rsd_num;

		if ( pose.glycan_tree_set() ) {
			prev_rsd_num = pose.glycan_tree_set()->get_parent( rsd_num );
		} else {
			prev_rsd_num = find_seqpos_of_saccharides_parent_residue( rsd );
		}

		Residue const & prev_rsd( pose.residue( prev_rsd_num ));
		// If this is not a saccharide residue, do nothing.
		if ( ! prev_rsd.is_carbohydrate() ) { return PREFERENCE_NA; }

		// If there is not an exocyclic bond between these two residues, this is not an omega about which we care.
		bool has_exocyclic_linkage;
		if ( pose.glycan_tree_set() ) {
			has_exocyclic_linkage = pose.glycan_tree_set()->has_exocyclic_glycosidic_linkage(rsd_num);
		} else {
			has_exocyclic_linkage = has_exocyclic_glycosidic_linkage(pose.conformation(), rsd_num );
		}


		if ( ! has_exocyclic_linkage ) { return PREFERENCE_NA; }

		// There is only a preference if the parent residue is an aldohexopyranose with an O4.
		chemical::carbohydrates::CarbohydrateInfoCOP info( prev_rsd.carbohydrate_info() );
		if ( info->is_aldose() && info->is_hexose() && info->is_pyranose() && prev_rsd.has( "O4" ) ) {
			// Get the atom indices we need.
			// We can assume that the ring we care about is the 1st ring, since this is a sugar.
			utility::vector1< uint > const ring_atoms( prev_rsd.type().ring_atoms( 1 ) );
			uint const O4( prev_rsd.atom_index( "O4" ) );

			switch ( is_atom_axial_or_equatorial_to_ring( prev_rsd, O4, ring_atoms ) ) {
			case AXIAL :
				return ANTI;
			case EQUATORIAL :
				return GAUCHE_EFFECT;
			case NEITHER :
				return ANTI;  // This is an assumption on my part.  ~Labonte
			}
		}
	}
	return PREFERENCE_NA;
}

/// @return  CHIEnergyFunctionLinkageType, which corresponds to one of six particular CHI Energy Functions specific to
/// the given glycosidic linkage, or LINKAGE_NA if the CHI Energy Function is not applicable to this torsion.
CHIEnergyFunctionLinkageType
get_CHI_energy_function_linkage_type_for_residue_in_pose(
	id::MainchainTorsionType torsion,
	pose::Pose const & pose,
	core::uint rsd_num )
{
	if ( torsion == id::phi_dihedral ) {
		return get_CHI_energy_function_linkage_type_for_phi_for_residue_in_pose( pose, rsd_num );
	} else if ( torsion == id::psi_dihedral ) {
		return get_CHI_energy_function_linkage_type_for_psi_for_residue_in_pose( pose, rsd_num );
	}
	return LINKAGE_NA;  // CHI Energy Functions only exist for phi and psi.
}

/// @author  Jared Adolf-Bryfogle <jadolfbr@gmail.com>
utility::vector1< CHIEnergyFunctionLinkageType >
get_linkage_types_for_dihedral( core::Size torsion )
{
	utility::vector1< CHIEnergyFunctionLinkageType > linkages( 2 );
	if ( torsion == id::phi_dihedral ) {
		linkages[ 1 ] = ALPHA_LINKS;
		linkages[ 2 ] = BETA_LINKS;
	} else if ( torsion == id::psi_dihedral ) {
		linkages[ 1 ] = _2AX_3EQ_4AX_LINKS;
		linkages[ 2 ] = _2EQ_3AX_4EQ_LINKS;

		//CHI in an an exocyclic 1-6 linkage.
		linkages.push_back( ALPHA6_LINKS );
		linkages.push_back( BETA6_LINKS );
	} else {
		//Should probably throw instead of exit.
		utility_exit_with_message( "no data for torsion" );
	}
	return linkages;
}

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core
