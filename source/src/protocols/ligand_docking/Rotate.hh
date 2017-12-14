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

#ifndef INCLUDED_protocols_ligand_docking_Rotate_hh
#define INCLUDED_protocols_ligand_docking_Rotate_hh

// Unit Headers
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/ligand_docking/DistributionMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//// Scoring grid headers

//// Project Headers
#include <core/kinematics/Jump.hh>
#include <core/grid/CartGrid.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

struct Rotate_info{ // including default values
	std::string chain;
	core::Size chain_id;// looking this up from chain is slow so we store it
	core::Size jump_id; // looking this up from chain is slow so we store it
	Distribution distribution;
	core::Size degrees;
	core::Size cycles;
	utility::vector1<core::Size> tag_along_chains; // must be one residue per chain, eg water, metal
	utility::vector1<core::Size> tag_along_jumps;
	utility::vector1<core::Size> tag_along_residues;
};

struct Ligand_info{
	core::conformation::ResidueCOPs residues; // multiple per main ligand
	core::conformation::ResidueCOPs tag_along_residues; // one per tag_along_chain
	int atr;
	int rep;
	core::kinematics::Jump jump;
	Ligand_info();
	Ligand_info(core::conformation::ResidueCOPs const & residues, int atr, int rep);
	Ligand_info(core::conformation::ResidueCOPs const & residues, std::pair<int,int> scores, core::kinematics::Jump jump);
	bool operator<(Ligand_info const & ligand_info) const;
	bool operator<(std::pair<int,int> const & scores) const;
	core::conformation::ResidueCOPs const & get_residues() const;
};

class Rotate: public protocols::moves::Mover
{
public:
	Rotate();
	Rotate(Rotate_info const & rotate_info); // moves Rotate_info, so by-value
	~Rotate() override;
	Rotate(Rotate const & that);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void rotate_ligand(
		utility::pointer::shared_ptr<core::grid::CartGrid<int> >  const & grid,
		core::pose::Pose & pose
	);

	void rotate_ligand(core::pose::Pose & pose);

	/// @brief  These should have repulsive and attractive scores under the threshold
	utility::vector1< Ligand_info> create_random_rotations(
		utility::pointer::shared_ptr<core::grid::CartGrid<int> > const & grid,
		protocols::rigid::RigidBodyMoverOP const mover,
		core::pose::Pose & pose,
		core::Size chain_begin
	)const;

	// utility::vector1<Ligand_info> create_random_rotations(
	//   protocols::rigid::RigidBodyMoverOP const ,
	//   core::Size const begin,
	//   core::pose::Pose & pose) const;

	Ligand_info create_random_rotation(
		utility::pointer::shared_ptr<core::grid::CartGrid<int> > const & grid,
		protocols::rigid::RigidBodyMoverOP const mover,
		core::Vector const & center,
		core::Size const begin,
		core::Size const end,
		core::pose::Pose & local_pose
	) const;
	/*
	Ligand_info create_random_rotation(
	protocols::rigid::RigidBodyMoverOP const mover,
	core::Vector const center,
	core::Size const begin,
	core::Size const end,
	core::pose::Pose & local_pose) const;
	*/

private:
	protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype_;
	Rotate_info rotate_info_;
}; // class Rotate

/// Convenience Functions for use with Rotate

bool check_score(
	Ligand_info const & ligand,
	core::Size const heavy_atom_number
);

bool check_RMSD(
	Ligand_info const & ligand,
	core::Size const heavy_atom_number,
	utility::vector1< Ligand_info> const & ligands
);

void apply_rotate(
	protocols::rigid::RigidBodyMoverOP mover,
	core::pose::Pose & pose,
	core::Vector const & center,
	core::Size jump_id,
	utility::vector1<core::Size> tag_along_chains
);

void add_ligand_conditionally(
	Ligand_info const & ligand_info,
	utility::vector1< Ligand_info> & ligands,
	core::Size const heavy_atom_number
);

} //namespace ligand_docking
} //namespace protocols

#endif
