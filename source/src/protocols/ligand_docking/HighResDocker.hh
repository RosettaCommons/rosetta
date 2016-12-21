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

#ifndef INCLUDED_protocols_ligand_docking_HighResDocker_hh
#define INCLUDED_protocols_ligand_docking_HighResDocker_hh

// Unit Headers
#include <protocols/ligand_docking/HighResDocker.fwd.hh>
#include <protocols/ligand_docking/ligand_options/Interface.fwd.hh>
#include <protocols/ligand_docking/MinimizeLigand.fwd.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <protocols/ligand_docking/TetherLigand.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/types.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers

// STL Headers
#include <map>
#include <string>
#include <vector>

#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class HighResDocker: public protocols::moves::Mover
{
public:

	HighResDocker();
	~HighResDocker() override;
	HighResDocker(HighResDocker const & that);
	HighResDocker(
		Size num_cycles,
		Size repack_every_Nth,
		core::scoring::ScoreFunctionOP score_fxn,
		MoveMapBuilderOP movemap_builder,
		std::string resfile=""
	);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	// XRW TEMP  std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;
	void apply(utility::vector1<core::pose::Pose> & poses, utility::vector1<core::Real> & current_scores, utility::vector1<char> qsar_chars, core::Size cycle);
	void apply(core::pose::Pose & pose, core::Real & current_score, char qsar_char, core::Size cycle);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Size num_cycles_;
	core::Size repack_every_Nth_;
	std::vector<std::string> chains_;
	core::scoring::ScoreFunctionOP score_fxn_;
	MoveMapBuilderOP movemap_builder_;
	std::string resfile_;

	MinimizeLigandOPs setup_ligands_to_minimize(core::pose::Pose pose, char chain = 0);
	TetherLigandOPs tether_ligands(core::pose::Pose & pose, char chain = 0);
	void remove_ligand_tethers(core::pose::Pose pose, TetherLigandOPs ligand_tethers);

	void enable_ligand_rotamer_packing(
		core::pose::Pose const & pose,
		core::Size const ligand_residue_id,
		core::pack::task::PackerTaskOP & pack_task
	) const;

	core::pack::task::PackerTaskOP
	make_packer_task(
		core::pose::Pose const & pose,
		bool all_residues= false
	) const;

	core::pack::task::PackerTaskOP
	make_packer_task_from_vector(
		core::pose::Pose const & pose,
		ligand_options::Interface const & allow_repack
	) const;

	utility::vector1<protocols::moves::MoverOP>
	create_rigid_body_movers(core::pose::Pose const & pose) const;

	void apply_rigid_body_moves(
		core::pose::Pose & pose,
		utility::vector1<protocols::moves::MoverOP> & rigid_body_movers
	);
};

//void setup_native_residue_favoring(core::pose::Pose & pose, core::pack::task::PackerTaskOP task);

} //namespace ligand_docking
} //namespace protocols

#endif
