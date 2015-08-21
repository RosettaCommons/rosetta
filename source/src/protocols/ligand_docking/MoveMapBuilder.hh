// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_MoveMapBuilder_hh
#define INCLUDED_protocols_ligand_docking_MoveMapBuilder_hh

// Unit Headers
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <protocols/ligand_docking/InterfaceBuilder.fwd.hh>

// Package Headers
#include <core/kinematics/MoveMap.fwd.hh>

//// Project Headers

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

/// @brief
class MoveMapBuilder: public utility::pointer::ReferenceCount
{
public:
	MoveMapBuilder();
	virtual ~MoveMapBuilder();
	MoveMapBuilder(MoveMapBuilder const & that);
	MoveMapBuilder(InterfaceBuilderOP sc, InterfaceBuilderOP bb, bool minimize_water);

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	core::kinematics::MoveMapOP
	build(core::pose::Pose const &) const;

	InterfaceBuilderOP
	get_sc_interface_builder()const;

	InterfaceBuilderOP
	get_bb_interface_builder()const;

	bool minimize_backbone();

private:
	InterfaceBuilderOP sc_interface_builder_; // which side chains to minimize?
	InterfaceBuilderOP bb_interface_builder_; // which backbone residues to minimize?

	bool minimize_water_;

	void
	set_all_chi(
		core::pose::Pose const & pose,
		core::kinematics::MoveMapOP movemap
	)const;

	void
	set_all_bb(
		core::pose::Pose const & pose,
		core::kinematics::MoveMapOP movemap
	)const;

};

} //namespace ligand_docking
} //namespace protocols

#endif
