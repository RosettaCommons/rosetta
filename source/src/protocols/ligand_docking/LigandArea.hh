// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/protocols/ligand_docking/ligand_options/Interface.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_LigandArea_hh
#define INCLUDED_protocols_ligand_docking_LigandArea_hh

//// Unit Headers
#include <protocols/ligand_docking/LigandArea.fwd.hh>

//// Package Headers
#include <core/types.hh>

//// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

//// C++ headers


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class LigandArea: public utility::pointer::ReferenceCount
{
public:
	virtual ~LigandArea();

	LigandArea();

	void parse_my_tag(
		utility::tag::TagCOP tag
	);

	char chain_;
	core::Real cutoff_;// angstroms from ligand to interface residue
	core::Real Calpha_restraints_;// size of one standard deviation (angstroms) for restraints on C-alphas
	core::Real minimize_ligand_; // size of one standard deviation (degrees) for ligand torsion angles
	core::Real tether_ligand_;
	core::Real high_res_angstroms_;
	core::Real high_res_degrees_;
	bool add_nbr_radius_;
	bool all_atom_mode_;
};

}
}

#endif
