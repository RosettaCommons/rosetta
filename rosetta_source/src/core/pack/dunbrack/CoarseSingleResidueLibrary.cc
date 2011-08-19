// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/CoarseSingleResidueLibrary.cc
/// @brief
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/pack/dunbrack/CoarseSingleResidueLibrary.hh>

// Package headers
#include <core/pack/dunbrack/CoarseRotamer.hh>


namespace core {
namespace pack {
namespace dunbrack {

CoarseSingleResidueLibrary::CoarseSingleResidueLibrary() {}

CoarseSingleResidueLibrary::~CoarseSingleResidueLibrary() {}

/// @brief Adheres to the contract from SingleLigandRotamerLibrary
Real
CoarseSingleResidueLibrary::rotamer_energy_deriv(
	conformation::Residue const & ,//rsd,
	RotamerLibraryScratchSpace & //scratch
) const
{
	return 0.0;
}

/// @brief Adheres to the contract from SingleLigandRotamerLibrary
Real
CoarseSingleResidueLibrary::rotamer_energy(
	conformation::Residue const & ,//rsd,
	RotamerLibraryScratchSpace & // scratch
) const
{
	return 0.0;
}

/// @brief Adheres to the contract from SingleLigandRotamerLibrary
Real
CoarseSingleResidueLibrary::best_rotamer_energy(
	conformation::Residue const & ,//rsd,
	bool ,//curr_rotamer_only,
	RotamerLibraryScratchSpace & // scratch
) const
{
	return 0.0;
}

/// @brief Adheres to the contract from SingleLigandRotamerLibrary
void
CoarseSingleResidueLibrary::fill_rotamer_vector(
	pose::Pose const &,// pose,
	scoring::ScoreFunction const &,// scorefxn,
	pack::task::PackerTask const &,// task,
	graph::GraphCOP,
	chemical::ResidueTypeCAP,// concrete_residue,
	conformation::Residue const&,// existing_residue,
	utility::vector1< utility::vector1< Real > > const &,// extra_chi_steps,
	bool,// buried,
	RotamerVector &// rotamers
) const
{}

/// @brief Adheres to the contract from SingleLigandRotamerLibrary
SingleResidueRotamerLibraryOP
CoarseSingleResidueLibrary::coarsify(
	coarse::Translator const & //map
) const
{
	return 0;
}

/// @brief Adheres to the contract from SingleLigandRotamerLibrary
void
CoarseSingleResidueLibrary::write_to_file(
	utility::io::ozstream & //out
) const
{}


void
CoarseSingleResidueLibrary::add_set(
	Size const,// phi_bin,
	Size const,// psi_bin,
	CoarseRotamerSetOP )//rotamers )
{}



} // namespace dunbrack
} // namespace scoring
} // namespace core

