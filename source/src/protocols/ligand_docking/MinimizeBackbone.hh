// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_MinimizeBackbone_hh
#define INCLUDED_protocols_ligand_docking_MinimizeBackbone_hh

// Unit Headers
#include <protocols/ligand_docking/MinimizeBackbone.fwd.hh>
#include <protocols/ligand_docking/ligand_options/Interface.fwd.hh>
#include <protocols/ligand_docking/InterfaceBuilder.fwd.hh>

// Package Headers
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>

//// Project Headers
#include <protocols/moves/Mover.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/id/AtomID.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

/// @brief
class MinimizeBackbone : public protocols::moves::Mover
{
public:
	MinimizeBackbone();
	MinimizeBackbone(InterfaceBuilderOP interface_builder);
	~MinimizeBackbone() override;
	MinimizeBackbone(MinimizeBackbone const & that);

	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	// map of ligand chains to minimize backbone around and how much minimization around each chain
	// Real is the Size of one standard deviation.  For restraints placed on C-alphas
	InterfaceBuilderOP interface_builder_;

	void reorder_foldtree_around_mobile_regions(
		ligand_options::Interface const & interface,
		core::pose::Pose & pose
	);

	core::kinematics::FoldTreeOP
	create_fold_tree_with_ligand_jumps_from_attach_pts(
		core::kinematics::FoldTree const &,
		ligand_options::Interface const & interface,
		core::pose::Pose & pose
	)const;

	core::kinematics::FoldTreeOP
	create_fold_tree_with_cutpoints(
		core::kinematics::FoldTreeCOP f,
		ligand_options::Interface const & interface,
		core::pose::Pose & pose
	);

	utility::vector1< protocols::loops::Loop >
	add_cut_points(
		core::kinematics::Edge const & edge,
		ligand_options::Interface const & interface,
		core::pose::Pose & pose
	);

	std::map<core::Size, core::Size> find_attach_pts(
		const ligand_options::Interface & interface,
		core::pose::Pose const & pose
	) const;

	void
	restrain_protein_Calphas(
		ligand_options::Interface const & interface,
		core::pose::Pose & pose
	);
	void restrain_protein_Calpha(
		ligand_options::Interface const & interface,
		core::pose::Pose & pose,
		core::Size residue_id,
		core::id::AtomID const & fixed_pt
	);
};

void restrict_to_protein_residues(
	ligand_options::Interface & interface,
	core::pose::Pose const & pose
);

void reorder_with_first_non_mobile_as_root(
	core::kinematics::FoldTreeOP f,
	const ligand_options::Interface & interface,
	core::pose::Pose & pose
);

core::Size find_attach_pt(
	core::Size const jump_id,
	ligand_options::Interface const & interface,
	core::pose::Pose const & pose
);

core::Size
find_peptide_attach_pt(
	int const & start,
	int const & stop,
	std::map<core::Size, core::Size > const & jump_to_attach
);


} //namespace ligand_docking
} //namespace protocols

#endif
