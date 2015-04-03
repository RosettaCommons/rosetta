// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/RandomConformerMover.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_RandomConformerMover_hh
#define INCLUDED_protocols_ligand_docking_RandomConformerMover_hh

#include <protocols/ligand_docking/RandomConformerMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {


/// @brief Replace the residue at the given position with a randomly selected
/// conformer from its rotamer library.
/// @details Only tested on ligand residues.
/// If using torsion restraints, should be wrapped in an UnconstrainedTorsionsMover.
class RandomConformerMover : public protocols::moves::Mover
{
public:

	RandomConformerMover(core::Size resid);
	virtual ~RandomConformerMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Size resid_;

}; // RandomConformerMover


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_RandomConformerMover_HH
