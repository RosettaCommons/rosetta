// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author James Thompson

#include <core/types.hh>

#include <core/id/NamedAtomID.hh>


#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/StealLigandMover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

StealLigandMover::StealLigandMover(
	core::pose::Pose const & source,
	core::id::NamedAtomID const & anchor_atom_dest,
	core::id::NamedAtomID const & anchor_atom_source,
	utility::vector1< core::id::NamedAtomID > const & ligand_indices
) :
	source_( source ),
	anchor_atom_dest_( anchor_atom_dest ),
	anchor_atom_source_( anchor_atom_source ),
	ligand_indices_( ligand_indices )
{}

void StealLigandMover::apply( core::pose::Pose & pose ) {
	steal_ligands(
		pose,
		source_,
		anchor_atom_dest_,
		anchor_atom_source_,
		ligand_indices_
	);
} // apply

std::string
StealLigandMover::get_name() const {
	return "StealLigandMover";
}


} // comparative_modeling
} // protocols
