// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/CoarseSingleResidueLibrary.hh
/// @brief
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_CoarseSingleResidueLibrary_hh
#define INCLUDED_core_pack_dunbrack_CoarseSingleResidueLibrary_hh

// Unit headers
#include <core/pack/dunbrack/CoarseSingleResidueLibrary.fwd.hh>

// Package headers
#include <core/pack/dunbrack/CoarseRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace pack {
namespace dunbrack {


class CoarseSingleResidueLibrary : public SingleResidueRotamerLibrary
{
public:

	CoarseSingleResidueLibrary();

	virtual ~CoarseSingleResidueLibrary();

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	Real
	rotamer_energy_deriv(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	Real
	rotamer_energy(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	Real
	best_rotamer_energy(
		conformation::Residue const & rsd,
		bool curr_rotamer_only,
		RotamerLibraryScratchSpace & scratch
	) const;

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	void
	fill_rotamer_vector(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCAP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		RotamerVector & rotamers
	) const;

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	SingleResidueRotamerLibraryOP
	coarsify(coarse::Translator const &map) const;

	/// @brief Adheres to the contract from SingleLigandRotamerLibrary
	virtual
	void
	write_to_file( utility::io::ozstream &out ) const;


	void
	add_set( Size const phi_bin, Size const psi_bin, CoarseRotamerSetOP rotamers );

private:

	ObjexxFCL::FArray2D< CoarseRotamerSetOP > rotamers_;

}; // CoarseSingleResidueLibrary


} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_CoarseSingleResidueLibrary_HH
