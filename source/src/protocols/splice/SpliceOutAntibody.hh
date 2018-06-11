// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/Splice.hh
/// @author Gideon Lapidoth (glapidoth@gmail.com)

#ifndef INCLUDED_protocols_splice_SpliceOutAntibody_hh
#define INCLUDED_protocols_splice_SpliceOutAntibody_hh

#include <protocols/splice/SpliceOutAntibody.fwd.hh>
#include <protocols/splice/SpliceSegment.fwd.hh>
#include <protocols/splice/SpliceManager.hh>
#include <protocols/splice/SpliceOutTail.hh>
#include <protocols/splice/Splice.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <protocols/task_operations/SeqprofConsensusOperation.fwd.hh>
#include <basic/database/open.hh>
#include <core/pack/task/PackerTask.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <protocols/moves/Mover.hh>


namespace protocols {
namespace splice {



class SpliceOutAntibody : public protocols::splice::SpliceOutTail
{

public:
	SpliceOutAntibody();
	void apply( Pose & pose ) override;
	virtual std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new SpliceOutAntibody ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )override;
	virtual ~SpliceOutAntibody();
	void find_disulfide_postions(core::pose::Pose const & pose, utility::vector1<core::Size> & cys_pos);
	void antibody_DB(std::string const s){antibody_DB_=s;}
	std::string antibody_DB(){return antibody_DB_;}
	void assign_from_res_to_res(core::pose::Pose const pose);
	void update_vl_vh_cut();
	void set_fold_tree_nodes(core::pose::Pose const & pose) override;
	void vl_vh_cut(core::Size i){vl_vh_cut_=i;}
	core::Size vl_vh_cut(){return vl_vh_cut_;}
	void handle_tail_mover_tag(TagCOP const tag,protocols::moves::Movers_map const & movers);
	void find_vl_vh_cut(core::pose::Pose pose);
	void set_source_from_to_res() override;
	void place_cut_site_in_segment(core::pose::Pose const & pose) override;
	core::Size set_anchor_res() override;
	void set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc) override;
	void superimpose_source_on_pose( core::pose::Pose const & pose, core::pose::Pose & source_pose ) override;
	void write_database_to_file(core::pose::Pose const & pose) override;
	void adjust_n_ter_tail_length(core::pose::Pose & pose);
	void build_ideal_segment(core::pose::Pose & pose) override;
	std::string name_for_filter() override;


	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string mover_name();


private:
	utility::vector1<core::Size> pose_cys_pos_, template_cys_pos_;//store pose and template cysteine positions. I assume that we are working with ScFv, VL before VH is sequence order.
	std::string antibody_DB_;
	core::Size vl_vh_cut_;
	protocols::moves::MoverOP tail_submover_;
	void (SpliceOut::*call_mover_tmp)(core::pose::Pose & pose,core::kinematics::MoveMapOP mm);

};

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_Splice_hh
