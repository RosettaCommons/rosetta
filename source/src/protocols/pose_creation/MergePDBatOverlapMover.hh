// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/MergePDBatOverlapMover.hh
/// @brief This class combines two pdbs with a known overlap
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_pose_creation_MergePDBatOverlapMover_hh
#define INCLUDED_protocols_pose_creation_MergePDBatOverlapMover_hh


#include <protocols/pose_creation/MergePDBatOverlapMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>


namespace protocols {
namespace pose_creation {

class MergePDBatOverlapMover : public moves::Mover {

public:

	MergePDBatOverlapMover();
	MergePDBatOverlapMover(core::scoring::ScoreFunctionOP sfxn);
	void increase_range_to_ignore_ss_element(core::pose::Pose const & pose, core::Size init_start, core::Size init_end, core::Size & ss_start, core::Size & ss_end);
	core::Size closest_non_overlap_residue(core::pose::Pose const & pose, core::Size resid, core::Size start_overlap_resid, core::Size end_overlap_resid);
	void merge_junction_sequence(Pose & pose,std::string pose_junction_seq,std::string attach_pose_junction_seq,core::Size first_overlap_position);
	bool merge_poses(Pose & pose,Pose & attach_pose);
	void assign_seq(Pose & pose, char residue_type, core::Size position);
	void minimize_overlap(Pose & pose,core::Size overlap_start,core::Size overlap_end);
	bool makeJunctions_apply(core::pose::Pose & pose, core::pose::Pose const & attach_pose, core::Size overlap_length, core::Real max_overlap_rmsd, std::string const & attachment_termini,std::string const & chain);
	bool apply_helper( core::pose::Pose & pose );
	void apply( core::pose::Pose & pose ) override;


	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	std::string get_name() const override;
	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
private:
	core::pose::PoseOP attach_pose_;
	core::Size overlap_length_;
	core::Real max_overlap_rmsd_;
	bool minimize_after_overlap_;
	std::string attachment_termini_;
	core::scoring::ScoreFunctionOP sfxn_;
	std::string attachment_chain_;
};

} // pose_creation
} // protocols

#endif

