// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ligand_docking/ProtLigEnsemble.hh
/// @brief  mover for docking ligand ensemble into protein ensemble
/// @author Darwin Fu

#ifndef INCLUDED_protocols_ligand_docking_ProtLigEnsemble_hh
#define INCLUDED_protocols_ligand_docking_ProtLigEnsemble_hh

// Unit Headers
#include <protocols/ligand_docking/ProtLigEnsemble.fwd.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/FinalMinimizer.hh>
#include <protocols/ligand_docking/ligand_options/Interface.fwd.hh>
#include <protocols/ligand_docking/MinimizeLigand.fwd.hh>
// AUTO-REMOVED #include <protocols/ligand_docking/MinimizeBackbone.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <protocols/ligand_docking/TetherLigand.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

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
#include <ObjexxFCL/FArray2D.hh>

// STL Headers
#include <map>
#include <string>
// AUTO-REMOVED #include <set>
#include <vector>

#include <utility/vector1.hh>
#include <utility/string_util.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

struct ProtLigPair_info{ // Store info for each protein-ligand pair

public:
	//change first three into vectors (create a vector of chains)
	core::Size mut_resid;
	char mut_target;
	char lig_chain;
	core::Real bind_data;
	bool has_bind;
	bool wild_type;

	std::string print_mut_string(){return (utility::to_string(mut_resid)+ "_" + utility::to_string(mut_target));};
	std::string print_string(){
		if ( wild_type ) { return ("WT_" + utility::to_string(lig_chain));}
		return (print_mut_string() + "_" + utility::to_string(lig_chain));
	};

	ProtLigPair_info(): mut_resid(0), mut_target(0), lig_chain(0), bind_data(0), has_bind(false), wild_type(false){};
	ProtLigPair_info(ProtLigPair_info const & ) = default;
};


class ProtLigEnsemble: public protocols::moves::Mover
{
public:

	ProtLigEnsemble();
	~ProtLigEnsemble() override;
	ProtLigEnsemble(ProtLigEnsemble const & that);

	virtual bool reinitialize_for_each_job() const override;
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

	// void enable_ligand_rotamer_packing(
	//  core::pose::Pose const & pose,
	//  core::Size const ligand_residue_id,
	//  core::pack::task::PackerTaskOP & pack_task
	// ) const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	utility::vector1<ProtLigPair_info> pose_infos_;
	core::Size ignore_correlation_until_;
	core::Size ignore_correlation_after_;
	core::Real distance_;
	core::Size num_cycles_;
	core::Size repack_every_Nth_;
	core::scoring::ScoreFunctionOP score_fxn_;


	core::scoring::ScoreFunctionOP final_score_fxn_;

	core::Real correlation_weight_;    //Weight term to adjust score by based on correlation

	utility::vector1<std::pair<core::Size, core::Real> > exp_ranks_; // Vector of chain ID/Affinity pairs

	utility::vector1<std::pair<core::Size, core::Real> > rosetta_lowest_scores_; // Vector of chain ID/Affinity pairs
	utility::vector1<core::pose::Pose> rosetta_lowest_poses_; // Vector of Poses


	void enable_ligand_rotamer_packing(
		core::pose::Pose const & pose,
		core::Size const ligand_residue_id,
		core::pack::task::PackerTaskOP & pack_task
	) const;

	void read_qsar_file(std::string filename);

	void prepare_single_ligand_pose(core::pose::Pose & pose, ProtLigPair_info & info, core::Size count);

	core::pack::task::PackerTaskOP make_packer_task(core::pose::Pose const & pose, utility::vector1<bool> const & interface_residues);

	core::Real qsar_correlation();
};

//non-member functions

ProtLigPair_info process_line(std::string & line);

bool sort_by_binding(
	const ProtLigPair_info & left_lig,
	const ProtLigPair_info & right_lig
);

void info_to_rank(utility::vector1<ProtLigPair_info> & vector);


} //namespace ligand_docking
} //namespace protocols

#endif
