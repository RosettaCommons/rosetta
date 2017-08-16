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


#include <protocols/simple_moves/MergePDBMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/pack/task/TaskFactory.hh>

namespace protocols {
namespace simple_moves {

class MergePDBMover : public moves::Mover {

public:
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
	MergePDBMover();
	utility::vector1<Overlap> determine_overlap(Pose const pose);
	utility::vector1<core::pose::PoseOP> generate_overlaps(Pose & pose, utility::vector1<MergePDBMover::Overlap> overlaps);
	utility::vector1<core::pose::PoseOP> pack_and_minimize(utility::vector1<core::pose::PoseOP> poses, utility::vector1<MergePDBMover::Overlap> overlaps,core::Real baseline_score);
	core::pose::PoseOP get_additional_output() override;

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
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
	Size overlap_scan_range_;
	core::Real design_range_;
	core::Real packing_range_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::pack::task::TaskFactoryOP task_factory_;
	bool do_minimize_;

};

} // simple_moves
} // protocols

#endif

