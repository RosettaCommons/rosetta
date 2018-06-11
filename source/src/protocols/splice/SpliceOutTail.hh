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

#ifndef INCLUDED_protocols_splice_SpliceOutTail_hh
#define INCLUDED_protocols_splice_SpliceOutTail_hh

#include <protocols/splice/SpliceOut.hh>
#include <protocols/splice/SpliceOutTail.fwd.hh>
#include <protocols/splice/SpliceSegment.fwd.hh>
#include <protocols/splice/SpliceManager.hh>
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
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>



namespace protocols {
namespace splice {



class SpliceOutTail : virtual public protocols::splice::SpliceOut{

public:
	SpliceOutTail();
	// SpliceOutTail(SpliceManager sm,std::string source_pdb, core::scoring::ScoreFunctionOP scorefxn, core::Real rms_cutoff,core::Real rms_cutoff_loop,
	//   core::pack::task::TaskFactoryOP task_factory, core::pack::task::TaskFactoryOP design_task_factory,
	//   bool poly_ala, core::pose::PoseOP source_pose,  core::kinematics::FoldTreeOP saved_fold_tree, bool design, bool allow_all_aa,
	//   bool thread_original_sequence, bool rtmin, std::string dbase_file_name, protocols::filters::FilterOP splice_filter,
	//   core::Size source_from_res, core::Real profile_weight_away_from_interface, bool restrict_to_repacking_chain2, core::Real design_shell,
	//   core::Real repack_shell, bool deleteHairpin, core::Size delete_hairpinC, core::Size delete_hairpinN,
	//   void(SpliceOut::*callMover)(core::pose::Pose & pose,core::kinematics::MoveMapOP mm),
	//   protocols::moves::MoverOP submover, bool write_to_database);
	virtual std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new SpliceOutTail ); }
	virtual ~SpliceOutTail();
	virtual void set_source_from_to_res() override;
	core::Size set_anchor_res() override;
	virtual void set_loop_length_change(protocols::protein_interface_design::movers::LoopLengthChange & llc) override;
	virtual void set_fold_tree_nodes(core::pose::Pose const & pose) override;
	virtual void place_cut_site_in_segment(core::pose::Pose const & pose) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string mover_name();
	void abstract_parse_tag(TagCOP const tag) override;
	void write_database_to_file(core::pose::Pose const & pose) override;
	void superimpose_source_on_pose( core::pose::Pose const & pose, core::pose::Pose & source_pose ) override;
	virtual void build_ideal_segment(core::pose::Pose & pose) override;
	std::string name_for_filter() override;


private:


};

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_Splice_hh
