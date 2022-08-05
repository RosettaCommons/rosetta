// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/dunbrack/RotamerLibraryScratchSpace.hh
/// @brief  Declaration of scratch space class for Dunbrack rotamer library
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>

// Package headers
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

namespace core {
namespace pack {
namespace dunbrack {


/// @details All the fixedsizearrays are allocated and initialized to 0
RotamerLibraryScratchSpace::RotamerLibraryScratchSpace() :
	fa_dun_tot_( 0.0 ),
	fa_dun_rot_( 0.0 ),
	fa_dun_semi_( 0.0 ),
	fa_dun_dev_( 0.0 )
{
}

RotamerLibraryScratchSpace::~RotamerLibraryScratchSpace() = default;

void
RotamerLibraryScratchSpace::extract_torsion_deriv(
	core::id::TorsionID const & tor_id,
	core::conformation::Residue const &rsd,
	core::pose::Pose const &pose,
	rotamers::TorsionEnergy & tderiv
) const {

	if ( tor_id.type() == id::BB ) {
		core::Size const scratch_index( get_scratch_index( tor_id, rsd, pose ) ); //The index of the relevant torsion, as indexed in the scratch space.
		// Note that the scratch space indexes only the relevant backbone torsions.  So, for example, if a residue type had
		// rotamers that were dependent on mainchain torsion indices 1, 2, and 4, the scratch space would have indices 1, 2,
		// and 3 corresponding to mainchain indices 1, 2, and 4.

		if ( scratch_index != 0 && scratch_index  <= dunbrack::DUNBRACK_MAX_BBTOR ) {
			tderiv.tot  += dE_dbb_[ scratch_index ];
			tderiv.dev  += dE_dbb_dev_[ scratch_index ];
			tderiv.rot  += dE_dbb_rot_[ scratch_index ];
			tderiv.semi += dE_dbb_semi_[ scratch_index ];
		}

	} else if ( tor_id.type() == id::CHI && tor_id.torsion() <= dunbrack::DUNBRACK_MAX_SCTOR ) {
		tderiv.tot  += dE_dchi_[ tor_id.torsion() ];
		tderiv.dev  += dE_dchi_dev_[ tor_id.torsion() ];
		tderiv.semi += dE_dchi_semi_[ tor_id.torsion() ];
	}
}

/// @brief Given a mainchain torsion index and a ResidueType, get the index of the corresponding torsion in the
/// data stored in the Dunbrack scratch space.
/// @details For most residue types, this just returns torsion_index.  The index is only different in cases in which
/// a residue type has rotamers that depend on a subset of mainchain torsions.  For example, if a residue's rotamers
/// depended on mainchain torsions 2, 3, and 4, then the scratch indices 1, 2, and 3 would correspond to mainchain
/// torsions 2, 3, and 4, respectively.  This function returns 0 if torsion_index is a torsion on which rotamers do
/// not depend.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::Size
RotamerLibraryScratchSpace::get_scratch_index(
	core::id::TorsionID const &torid,
	core::conformation::Residue const &rsd,
	core::pose::Pose const &pose
) const {
	// Special-case logic for peptoids:
	if ( rsd.is_peptoid() ) {
		if ( torid.rsd() == rsd.connected_residue_at_lower() && torid.torsion() == pose.residue( rsd.connected_residue_at_lower() ).mainchain_torsions().size() ) {
			return 3;
		} else if ( torid.rsd() == rsd.seqpos() && torid.torsion() <=2 ) {
			return torid.torsion();
		} else {
			return 0;
		}
	}

	// In general, torsion id must be within this residue:
	if ( torid.rsd() != rsd.seqpos() ) return 0;

	core::Size scratch_index(torid.torsion());

	// Figure out whether this is a residue type with rotamers that depend on only a subset of mainchain torsion angles:
	core::chemical::rotamers::RotamerLibrarySpecificationCOP rotlibspec( rsd.type().rotamer_library_specification() );
	if ( rotlibspec != nullptr ) {
		core::chemical::rotamers::NCAARotamerLibrarySpecificationCOP ncaa_rotlibspec( utility::pointer::dynamic_pointer_cast< core::chemical::rotamers::NCAARotamerLibrarySpecification const>( rotlibspec ) );
		if ( ncaa_rotlibspec != nullptr ) {
			debug_assert( ncaa_rotlibspec->rotamer_bb_torsion_indices().size() < rsd.mainchain_torsions().size() ); //Should always be true.
			if ( ncaa_rotlibspec->rotamer_bb_torsion_indices().size() < rsd.mainchain_torsions().size() - 1 ) {
				scratch_index = 0;
				for ( core::Size i(1), imax(ncaa_rotlibspec->rotamer_bb_torsion_indices().size()); i<=imax; ++i ) {
					if ( ncaa_rotlibspec->rotamer_bb_torsion_indices()[i] == torid.torsion() ) scratch_index = i;
				}
			}
		}
	}

	return scratch_index;
}

}
}
}

