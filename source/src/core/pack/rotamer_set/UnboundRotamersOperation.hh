// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/UnboundRotamersOperation.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_core_pack_rotamer_set_UnboundRotamersOperation_hh
#define INCLUDED_core_pack_rotamer_set_UnboundRotamersOperation_hh

#include <core/pack/rotamer_set/UnboundRotamersOperation.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#ifdef WIN32
#include <core/pose/Pose.hh>
#else
#include <core/pose/Pose.fwd.hh>
#endif

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {


/// @brief Adds in rotamers from the "unbound" or native structure(s),
/// pulled from one or more PDBs supplied on the command line.
/// @details Sequence numbering matters -- rotamers will only be added
/// if sequence numbers match AND the ResidueType is allowed by the PackerTask.
/// By itself, this class does NOT grant a Dunbrack energy bonus to the native rotamer(s).
class UnboundRotamersOperation : public core::pack::rotamer_set::RotamerSetOperation
{
public:

	UnboundRotamersOperation();
	virtual ~UnboundRotamersOperation();

	/// @brief Adds rotamers from the specified pose to the unbound collection.
	virtual void add_pose(core::pose::PoseCOP pose);

	virtual Size total_residue();

	/// @brief Loads poses from the -unboundrot flag.
	virtual void initialize_from_command_line();

	virtual
	core::pack::rotamer_set::RotamerSetOperationOP
	clone() const;

	virtual
	void
	alter_rotamer_set(
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		task::PackerTask const & ptask,
		graph::GraphCOP packer_neighbor_graph,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);

private:
	Size total_rot_;
	utility::vector1< core::pose::PoseCOP > poses_;

}; // UnboundRotamersOperation


} // namespace rotamer_set
} // namespace pack
} // namespace core

#endif // INCLUDED_core_pack_rotamer_set_UnboundRotamersOperation_HH
