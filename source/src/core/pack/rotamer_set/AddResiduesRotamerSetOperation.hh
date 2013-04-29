// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file core/pack/rotamer_set/AddResiduesRotamersOperation.hh
/// @brief Transform a vector of residues onto the backbone of the appropriate pose position and
/// add the transformed residues to the rotamer set
/// @author Tim Jacobs

#ifndef INCLUDED_core_pack_rotamer_set_AddResiduesRotamerSetOperation_hh
#define INCLUDED_core_pack_rotamer_set_AddResiduesRotamerSetOperation_hh

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/Residue.hh>

namespace core {
namespace pack {
namespace rotamer_set {

///@brief Adds in rotamers from a list of Residues,
class AddResiduesRotamerSetOperation : public core::pack::rotamer_set::RotamerSetOperation
{
public:

AddResiduesRotamerSetOperation(
	utility::vector1<core::conformation::ResidueOP> const & residues
);

core::pack::rotamer_set::RotamerSetOperationOP
clone() const;

void
alter_rotamer_set(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & /*sfxn*/,
	core::pack::task::PackerTask const & /*ptask*/,
	core::graph::GraphCOP /*packer_neighbor_graph*/,
	core::pack::rotamer_set::RotamerSet & rotamer_set
);

private:
	utility::vector1<core::conformation::ResidueOP> residues_;
};

} //rotamer_set
} //pack
} //core

#endif /* core_pack_rotamer_set_AddResiduesRotamerSetOperation_hh */
