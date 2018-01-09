// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/MergePDBatOverlapMover.hh
/// @brief This class combines two pdbs with a known overlap
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_MergePDBatOverlapMover_hh
#define INCLUDED_protocols_simple_moves_MergePDBatOverlapMover_hh


#include <protocols/simple_moves/MergePDBatOverlapMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/pack/task/TaskFactory.hh>

namespace protocols {
namespace simple_moves {

class MergePDBatOverlapMover : public moves::Mover {

public:

	MergePDBatOverlapMover();
	MergePDBatOverlapMover(core::scoring::ScoreFunctionOP sfxn);
	void increase_range_to_ignore_ss_element(core::pose::Pose const & pose, Size init_start, Size init_end, Size & ss_start, Size & ss_end);
	Size closest_non_overlap_residue(core::pose::Pose const & pose, core::Size resid, core::Size start_overlap_resid, core::Size end_overlap_resid);
	void merge_junction_sequence(Pose & pose,std::string pose_junction_seq,std::string attach_pose_junction_seq,Size first_overlap_position);
	bool merge_poses(Pose & pose,Pose & attach_pose);
	void assign_seq(Pose & pose, char residue_type, Size position);
	void minimize_overlap(Pose & pose,Size overlap_start,Size overlap_end);
	bool makeJunctions_apply(core::pose::Pose & pose, core::pose::Pose const & attach_pose, Size overlap_length, core::Real max_overlap_rmsd, std::string attachment_termini);
	bool apply_helper( core::pose::Pose & pose );
	void apply( core::pose::Pose & pose ) override;


	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string get_name() const override;
	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
private:
	core::pose::PoseOP attach_pose_;
	Size overlap_length_;
	core::Real max_overlap_rmsd_;
	bool minimize_after_overlap_;
	std::string attachment_termini_;
	core::scoring::ScoreFunctionOP sfxn_;

};

} // simple_moves
} // protocols

#endif

