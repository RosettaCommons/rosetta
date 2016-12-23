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

#ifndef INCLUDED_protocols_ligand_docking_HighResEnsemble_hh
#define INCLUDED_protocols_ligand_docking_HighResEnsemble_hh

// Unit Headers
#include <protocols/ligand_docking/HighResEnsemble.fwd.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/FinalMinimizer.hh>
#include <protocols/ligand_docking/ligand_options/Interface.fwd.hh>
#include <protocols/ligand_docking/MinimizeLigand.fwd.hh>
// AUTO-REMOVED #include <protocols/ligand_docking/MinimizeBackbone.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <protocols/ligand_docking/TetherLigand.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <protocols/ligand_docking/ResidueTorsionRestraints.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/types.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
//#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers

// STL Headers
#include <map>
#include <string>
// AUTO-REMOVED #include <set>
#include <vector>

#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

//struct HighRes_info{
//
//public:
// //Basic info
// core::Size jump_id;
// core::Real exp_data;
//
// //Rosetta scores and poses
// core::pose::Pose old_pose; //Pose for the chain in Rosetta
// core::pose::Pose current_pose; //Pose for the chain in Rosetta
// core::pose::Pose lowest_pose; //Pose for the chain in Rosetta
//
//
//
// core::Real current_score; //Score for the chain in Rosetta
// core::Real old_score; //Score for the chain in Rosetta
// core::Real lowest_score; //Score for the chain in Rosetta
//
//};

class HighResEnsemble: public protocols::moves::Mover
{
public:

	HighResEnsemble();
	~HighResEnsemble() override;
	HighResEnsemble(HighResEnsemble const & that);

	virtual bool reinitialize_for_each_job() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP const tag,
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
	core::Size num_cycles_;
	core::Size repack_every_Nth_;
	std::vector<std::string> chains_;
	core::scoring::ScoreFunctionOP score_fxn_;
	MoveMapBuilderOP movemap_builder_;
	std::string resfile_;

	//variables for final minimize mover
	core::scoring::ScoreFunctionOP final_score_fxn_;
	MoveMapBuilderOP final_movemap_builder_;

	core::Real correlation_weight_;    //Weight term to adjust score by based on correlation

	utility::vector1<std::pair<core::Size, core::Real> > exp_ranks_; // Vector of chain ID/Affinity pairs

	utility::vector1<std::pair<core::Size, core::Real> > rosetta_current_scores_; // Vector of chain ID/Affinity pairs
	utility::vector1<std::pair<core::Size, core::Real> > rosetta_lowest_scores_; // Vector of chain ID/Affinity pairs
	utility::vector1<std::pair<core::Size, core::Real> > rosetta_old_scores_; // Vector of chain ID/Affinity pairs


	utility::vector1<core::pose::Pose> rosetta_current_poses_; // Vector of Poses
	utility::vector1<core::pose::Pose> rosetta_old_poses_; // Vector of Poses
	utility::vector1<core::pose::Pose> rosetta_lowest_poses_; // Vector of Poses


	utility::vector1<std::string> rosetta_names_;    //QSAR compound names from file only
	utility::vector1<char> rosetta_chars_;          //Ligand chain identifiers

	void prepare_single_ligand_pose(core::pose::Pose pose, core::Size chain_to_keep);

	// void get_ranks_from_file(core::pose::Pose const pose, std::string filename);

	core::Real qsar_correlation();
};

//non-member functions

bool sort_by_second(
	std::pair<core::Size, core::Real> left_lig,
	std::pair<core::Size, core::Real> right_lig
);

core::Real spearman(
	utility::vector1<std::pair<core::Size, core::Real> > vector_exp,
	utility::vector1<std::pair<core::Size, core::Real> > vector_rosetta
);
void vector_to_rank(utility::vector1<std::pair<core::Size, core::Real> > & vector);

//void setup_native_residue_favoring(core::pose::Pose & pose, core::pack::task::PackerTaskOP task);

} //namespace ligand_docking
} //namespace protocols

#endif
