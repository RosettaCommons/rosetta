// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/MergePDBMover.hh
/// @brief This class will allign & combine parts of the pdb.
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_MergePDBMover_hh
#define INCLUDED_protocols_simple_moves_MergePDBMover_hh

#include <basic/datacache/DataMap.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MergePDBMover.fwd.hh>
#include <utility/vector1.hh>

#include <core/pack/task/TaskFactory.hh>

namespace protocols {
namespace simple_moves {

class MergePDBMover : public moves::Mover {

public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;

	struct Overlap{
		Size start_overlap_xmlPose;
		Size end_overlap_xmlPose;
		Size start_overlap_pose;
		Size end_overlap_pose;
		Overlap(Size start_overlap_xmlPose_i, Size end_overlap_xmlPose_i, Size start_overlap_pose_i, Size end_overlap_pose_i){
			start_overlap_xmlPose = start_overlap_xmlPose_i;
			end_overlap_xmlPose = end_overlap_xmlPose_i;
			start_overlap_pose = start_overlap_pose_i;
			end_overlap_pose = end_overlap_pose_i;
		}
		void print(){
			std::cout << start_overlap_xmlPose <<"," << end_overlap_xmlPose << "::" << start_overlap_pose << "," << end_overlap_pose << std::endl;
		}
	};
	/// @brief  Constructor
	MergePDBMover();
	/// @brief Determines the overlaps. stores the start and end position in the struct Overlap
	utility::vector1<Overlap> determine_overlap(Pose const pose,Size chain_id);
	/// @brief fast check to maake sure the pose isn't a structural duplicate
	bool check_duplicate(Pose &pose, utility::vector1<core::pose::PoseOP> outputPoses);
	/// @brief Uses the overlap to generate poses
	utility::vector1<core::pose::PoseOP> generate_overlaps(Pose & pose, utility::vector1<MergePDBMover::Overlap> overlaps,Size chain_id);
	/// @brief Figures out the closest residue that's not part of the overlap.
	core::Size closest_non_overlap_residue(core::pose::Pose const & pose, core::Size resid, core::Size start_overlap_resid, core::Size end_overlap_resid);
	/// @brief Gets the entire SS element the match is on
	void increase_range_to_ignore_ss_element(core::pose::Pose const & pose, Size init_start, Size init_end, Size & ss_start, Size & ss_end);
	/// @brief Copies the sequence in the overlap region as appropriate. This sets the initial residues before the pack and minimize is called
	void copy_sequence(core::Size start_overlap_resid, core::Size end_overlap_resid, core::Size start_overlap_input_pose_resid,core::Size start_overlap_xml_pose_resid,core::pose::Pose const & input_pose,core::pose::Pose const & xml_pose,core::pose::Pose & output_pose);
	/// @brief packs and minimizes if no clashes as determined by score0
	utility::vector1<core::pose::PoseOP> pack_and_minimize(utility::vector1<core::pose::PoseOP> poses,core::Real baseline_score);
	/// @brief any poses that score lower than the input files is output
	core::pose::PoseOP get_additional_output() override;

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	utility::vector1<core::pose::PoseOP> outputPoses_;
	utility::vector1<bool> outputYet_;
	core::pose::PoseOP xml_input_pose_;
	std::string overlap_location_pose_;
	core::Real overlap_max_rmsd_;
	Size overlap_length_;
	Size overlap_scan_range_cmdLine_;
	Size overlap_scan_range_xml_;
	core::Real design_range_;
	core::Real packing_range_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::scoring::ScoreFunctionOP asymm_score_;
	core::pack::task::TaskFactoryOP task_factory_;
	bool do_minimize_;
	std::string chain_;
	std::string symm_file_;
	std::string no_design_label_;
	std::string init_overlap_sequence_;
	core::Real duplicate_rmsd_pose_threshold_;
	bool detect_disulf_before_repack_;
};

} // simple_moves
} // protocols

#endif

