// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/DunbrackEnergy.cc
/// @brief  Dunbrack energy method implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/dunbrack/DunbrackEnergy.hh>
#include <core/pack/dunbrack/DunbrackEnergyCreator.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/id/TorsionID.hh>

// Utility headers
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

using namespace scoring;
using namespace scoring::methods;

/// @details This must return a fresh instance of the DunbrackEnergy class,
/// never an instance already in use
scoring::methods::EnergyMethodOP
DunbrackEnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const &
) const {
	return scoring::methods::EnergyMethodOP( new DunbrackEnergy );
}

scoring::ScoreTypes
DunbrackEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_dun );
	sts.push_back( fa_dun_dev );
	sts.push_back( fa_dun_rot );
	sts.push_back( fa_dun_semi );
	return sts;
}


/// ctor
DunbrackEnergy::DunbrackEnergy() :
	parent( EnergyMethodCreatorOP( new DunbrackEnergyCreator ) )
{}

DunbrackEnergy::~DunbrackEnergy() = default;

/// clone
EnergyMethodOP
DunbrackEnergy::clone() const
{
	return EnergyMethodOP( new DunbrackEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

/// @details Allocate the scratch space object on the stack to
/// alieviate thread-safety concerns.  Scratch does not use new.
void
DunbrackEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) )  return;

	if ( rsd.is_virtual_residue() ) return;

	//Returns the equivalent L-amino acid library if a D-amino acid is provided
	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib = rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );

	if ( ! rotlib || rsd.has_variant_type( core::chemical::SC_BRANCH_POINT ) ) return;

	dunbrack::RotamerLibraryScratchSpace scratch;
	emap[ fa_dun ] += rotlib->rotamer_energy( rsd, scratch );
	emap[ fa_dun_rot ] += scratch.fa_dun_rot();
	emap[ fa_dun_semi ] += scratch.fa_dun_semi();
	emap[ fa_dun_dev  ] += scratch.fa_dun_dev();
}

bool DunbrackEnergy::defines_dof_derivatives( pose::Pose const & ) const { return true; }

Real
DunbrackEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & ,//min_data,
	id::DOF_ID const & ,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & ,//pose,
	scoring::ScoreFunction const & ,//sfxn,
	scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return 0.0;

	Real deriv( 0.0 );
	Real deriv_dev( 0.0 );
	Real deriv_rot( 0.0 );
	Real deriv_semi( 0.0 );
	if ( ! tor_id.valid() ) return 0.0;

	debug_assert( rsd.seqpos() == tor_id.rsd() );

	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib =
		rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );
	if ( ! rsd.is_protein() || ! rotlib || rsd.has_variant_type( core::chemical::SC_BRANCH_POINT ) )  { return 0.0; }

	dunbrack::RotamerLibraryScratchSpace scratch;
	rotlib->rotamer_energy_deriv( rsd, scratch );
	if ( tor_id.type() == id::BB ) {
		core::Size const scratch_index( get_scratch_index( tor_id.torsion(), rsd ) ); //The index of the relevant torsion, as indexed in the scratch space.
		// Note that the scratch space indexes only the relevant backbone torsions.  So, for example, if a residue type had
		// rotamers that were dependent on mainchain torsion indices 1, 2, and 4, the scratch space would have indices 1, 2,
		// and 3 corresponding to mainchain indices 1, 2, and 4.

		if ( scratch_index == 0 ) {
			deriv = 0;
			deriv_dev = 0;
			deriv_rot = 0;
			deriv_semi = 0;
		} else if ( scratch_index  <= DUNBRACK_MAX_BBTOR ) {
			deriv      = scratch.dE_dbb()[ scratch_index ];
			deriv_dev  = scratch.dE_dbb_dev()[ scratch_index ];
			deriv_rot  = scratch.dE_dbb_rot()[ scratch_index ];
			deriv_semi = scratch.dE_dbb_semi()[ scratch_index ];
		}
	} else if ( tor_id.type() == id::CHI && tor_id.torsion() <= dunbrack::DUNBRACK_MAX_SCTOR ) {
		deriv      = scratch.dE_dchi()[ tor_id.torsion() ];
		deriv_dev  = scratch.dE_dchi_dev()[ tor_id.torsion() ];
		deriv_semi = scratch.dE_dchi_semi()[ tor_id.torsion() ];
	}

	return numeric::conversions::degrees( weights[ fa_dun ] * deriv + weights[ fa_dun_dev ] * deriv_dev + weights[ fa_dun_rot ] * deriv_rot + weights[ fa_dun_semi ] * deriv_semi);
}


Real
DunbrackEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	scoring::ScoreFunction const &,
	scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	Real deriv( 0.0 );
	Real deriv_dev( 0.0 );
	Real deriv_rot( 0.0 );
	Real deriv_semi( 0.0 );
	if ( !tor_id.valid() )  return 0.0;

	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib =
		rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue( tor_id.rsd() ).type() );

	if ( pose.residue( tor_id.rsd() ).is_virtual_residue() ) return 0.0;

	/// ASSUMPTION: Derivatives for amino acids only!
	if ( ! rotlib || ! pose.residue_type( tor_id.rsd() ).is_protein() || pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::SC_BRANCH_POINT) )  { return 0.0; }

	dunbrack::RotamerLibraryScratchSpace scratch;
	rotlib->rotamer_energy_deriv( pose.residue( tor_id.rsd() ), scratch );

	if ( tor_id.type() == id::BB ) {
		core::Size const scratch_index( get_scratch_index( tor_id.torsion(), pose.residue(tor_id.rsd()) ) ); //The index of the relevant torsion, as indexed in the scratch space.
		// Note that the scratch space indexes only the relevant backbone torsions.  So, for example, if a residue type had
		// rotamers that were dependent on mainchain torsion indices 1, 2, and 4, the scratch space would have indices 1, 2,
		// and 3 corresponding to mainchain indices 1, 2, and 4.

		if ( scratch_index == 0 ) {
			deriv = 0;
			deriv_dev = 0;
			deriv_rot = 0;
			deriv_semi = 0;
		} else if ( scratch_index  <= dunbrack::DUNBRACK_MAX_BBTOR ) {
			deriv      = scratch.dE_dbb()[ scratch_index ];
			deriv_dev  = scratch.dE_dbb_dev()[ scratch_index ];
			deriv_rot  = scratch.dE_dbb_rot()[ scratch_index ];
			deriv_semi = scratch.dE_dbb_semi()[ scratch_index ];
		}
	} else if ( tor_id.type() == id::CHI && tor_id.torsion() <= dunbrack::DUNBRACK_MAX_SCTOR ) {
		deriv      = scratch.dE_dchi()[ tor_id.torsion() ];
		deriv_dev  = scratch.dE_dchi_dev()[ tor_id.torsion() ];
		deriv_semi = scratch.dE_dchi_semi()[ tor_id.torsion() ];
	}

	return numeric::conversions::degrees( weights[ fa_dun ] * deriv + weights[ fa_dun_dev ] * deriv_dev + weights[ fa_dun_rot ] * deriv_rot + weights[ fa_dun_semi ] * deriv_semi );
}

/// @brief DunbrackEnergy is context independent; indicates that no context graphs are required
void
DunbrackEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

/// @brief Given a mainchain torsion index and a ResidueType, get the index of the corresponding torsion in the
/// data stored in the Dunbrack scratch space.
/// @details For most residue types, this just returns torsion_index.  The index is only different in cases in which
/// a residue type has rotamers that depend on a subset of mainchain torsions.  For example, if a residue's rotamers
/// depended on mainchain torsions 2, 3, and 4, then the scratch indices 1, 2, and 3 would correspond to mainchain
/// torsions 2, 3, and 4, respectively.  This function returns 0 if torsion_index is a torsion on which rotamers do
/// not depend.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::Size
DunbrackEnergy::get_scratch_index(
	core::Size const torsion_index,
	core::conformation::Residue const &rsd
) const {
	core::Size scratch_index(torsion_index);

	// Figure out whether this is a residue type with rotamers that depend on only a subset of mainchain torsion angles:
	core::chemical::rotamers::RotamerLibrarySpecificationCOP rotlibspec( rsd.type().rotamer_library_specification() );
	if ( rotlibspec != nullptr ) {
		core::chemical::rotamers::NCAARotamerLibrarySpecificationCOP ncaa_rotlibspec( utility::pointer::dynamic_pointer_cast< core::chemical::rotamers::NCAARotamerLibrarySpecification const>( rotlibspec ) );
		if ( ncaa_rotlibspec != nullptr ) {
			debug_assert( ncaa_rotlibspec->rotamer_bb_torsion_indices().size() < rsd.mainchain_torsions().size() ); //Should always be true.
			if ( ncaa_rotlibspec->rotamer_bb_torsion_indices().size() < rsd.mainchain_torsions().size() - 1 ) {
				scratch_index = 0;
				for ( core::Size i(1), imax(ncaa_rotlibspec->rotamer_bb_torsion_indices().size()); i<=imax; ++i ) {
					if ( ncaa_rotlibspec->rotamer_bb_torsion_indices()[i] == torsion_index ) scratch_index = i;
				}
			}
		}
	}

	return scratch_index;
}

core::Size
DunbrackEnergy::version() const
{
	return 1; // Initial versioning
}


} // dunbrack
} // pack
} // core

