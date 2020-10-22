// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/dunbrack/ResidueDOFReporter.hh
/// @brief   Class to measure the DOFs used by a RotamerLibrary
/// @author  Andrew Leaver-Fay

// Unit Headers
#include <core/pack/dunbrack/StandardDOFReporters.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>


namespace core {
namespace pack {
namespace dunbrack {

MainchainTorsionReporter::MainchainTorsionReporter(
	Size tor_ind,
	Size upper_tor_ind,
	Real neutral_val,
	bool flip_neutral_for_d_aa,
	bool flip_neutral_for_mirrored
) :
	tor_ind_( tor_ind ),
	upper_tor_ind_( upper_tor_ind ),
	neutral_val_( neutral_val ),
	flip_neutral_for_d_aa_( flip_neutral_for_d_aa ),
	flip_neutral_for_mirrored_( flip_neutral_for_mirrored )
{}


MainchainTorsionReporter::~MainchainTorsionReporter() = default;

Real
MainchainTorsionReporter::get_dof(
	conformation::Residue const & rsd,
	pose::Pose const &
) const
{
	Real const d_multiplier = ( flip_neutral_for_d_aa_ ?
		(rsd.type().is_d_aa() ? -1.0 : 1.0) :
		(flip_neutral_for_mirrored_ ?
		(rsd.type().is_mirrored_type() ? -1.0 : 1.0) :
		1.0));
	if ( ( rsd.is_lower_terminus() || !rsd.has_lower_connect() ) && tor_ind_ == 1 ) {
		return d_multiplier * neutral_val_;
	} else if ( ( rsd.is_upper_terminus() || !rsd.has_upper_connect() ) && tor_ind_ == upper_tor_ind_ ) {
		return d_multiplier * neutral_val_;
	} else if ( tor_ind_ == 0 ) {
		return 0.0;
	}

	return rsd.mainchain_torsion( tor_ind_ );
}

void
MainchainTorsionReporter::insert_atoms_defining_dof(
	conformation::Residue const & rsd,
	pose::Pose const &,
	std::set< id::PartialAtomID > & atoms
) const
{
	if (
			! ( (rsd.is_lower_terminus() || !rsd.has_lower_connect()) && tor_ind_ == 1 ) &&
			! ( (rsd.is_upper_terminus() || !rsd.has_upper_connect()) && tor_ind_ == upper_tor_ind_ ) ) {
		conformation::insert_partial_atom_ids_for_mainchain_torsion(
			rsd, tor_ind_, atoms );
	}
}


PeptideTorsionReporter::PeptideTorsionReporter(
	Size tor_ind,
	Size upper_tor_ind,
	Real neutral_val,
	bool flip_neutral_for_mirrored
) :
	tor_ind_( tor_ind ),
	upper_tor_ind_( upper_tor_ind ),
	neutral_val_( neutral_val ),
	flip_neutral_for_mirrored_( flip_neutral_for_mirrored )
{}


PeptideTorsionReporter::~PeptideTorsionReporter() = default;

Real
PeptideTorsionReporter::get_dof(
	conformation::Residue const & rsd,
	pose::Pose const &
) const
{
	int bb_torsion_index = bb_torsion_index_for_rsd( rsd );

	Real const d_multiplier = ( flip_neutral_for_mirrored_ ?
		(rsd.type().is_mirrored_type() ? -1.0 : 1.0) : 1.0);
	if ( ( rsd.is_lower_terminus() || !rsd.has_lower_connect() ) && tor_ind_ == 1 ) {
		return d_multiplier * neutral_val_;
	} else if ( ( rsd.is_upper_terminus() || !rsd.has_upper_connect() ) && tor_ind_ == upper_tor_ind_ ) {
		return d_multiplier * neutral_val_;
	}
	if ( bb_torsion_index != 0 ) {
		return rsd.mainchain_torsion( bb_torsion_index );
	}
	return 0.0;
}

void
PeptideTorsionReporter::insert_atoms_defining_dof(
	conformation::Residue const & rsd,
	pose::Pose const &,
	std::set< id::PartialAtomID > & atoms
) const
{
	int bb_torsion_index = bb_torsion_index_for_rsd( rsd );
	if ( ( rsd.is_lower_terminus() || !rsd.has_lower_connect() ) && tor_ind_ == 1 ) {
		return;
	} else if ( ( rsd.is_upper_terminus() || !rsd.has_upper_connect() ) && tor_ind_ == upper_tor_ind_ ) {
		return;
	}
	if ( bb_torsion_index != 0 ) {
		conformation::insert_partial_atom_ids_for_mainchain_torsion(
			rsd, bb_torsion_index, atoms );
	}
}

int
PeptideTorsionReporter::bb_torsion_index_for_rsd( conformation::Residue const & rsd ) const
{
	int bb_torsion_index(-1);
	// Figure out whether this is a residue type with rotamers that depend on only a subset of mainchain torsion angles:
	core::chemical::rotamers::RotamerLibrarySpecificationCOP rotlibspec( rsd.type().rotamer_library_specification() );
	if ( rotlibspec != nullptr ) {
		core::chemical::rotamers::NCAARotamerLibrarySpecificationCOP ncaa_rotlibspec( utility::pointer::dynamic_pointer_cast< core::chemical::rotamers::NCAARotamerLibrarySpecification const>( rotlibspec ) );
		if ( ncaa_rotlibspec != nullptr ) {
			if ( tor_ind_ <= ncaa_rotlibspec->rotamer_bb_torsion_indices().size() ) {
				bb_torsion_index = static_cast< int >( ncaa_rotlibspec->rotamer_bb_torsion_indices()[tor_ind_] );
			} else {
				bb_torsion_index = 0;
			}
		}
	}
	if ( bb_torsion_index == -1 ) { //If no rotlibspec, assume all mainchain torsions in this residue are relevant.
		bb_torsion_index = tor_ind_;
	}
	return bb_torsion_index;
}

}
}
}

