// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_StealLigandMover_hh
#define INCLUDED_protocols_comparative_modeling_StealLigandMover_hh

// Unit headers
#include <protocols/comparative_modeling/StealLigandMover.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/NamedAtomID.hh>

#include <protocols/moves/Mover.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

class StealLigandMover : public protocols::moves::Mover {

public:
	StealLigandMover(
		core::pose::Pose const & source,
		core::id::NamedAtomID const & anchor_atom_dest,
		core::id::NamedAtomID const & anchor_atom_source,
		utility::vector1< core::id::NamedAtomID > const & ligand_indices
	);

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	core::pose::Pose const & source_;
	core::id::NamedAtomID const anchor_atom_dest_;
	core::id::NamedAtomID const anchor_atom_source_;
	utility::vector1< core::id::NamedAtomID > const ligand_indices_;
};

} // comparative_modeling
} // protocols

#endif
