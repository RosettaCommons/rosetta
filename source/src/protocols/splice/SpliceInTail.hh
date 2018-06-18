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

#ifndef INCLUDED_protocols_splice_SpliceInTail_hh
#define INCLUDED_protocols_splice_SpliceInTail_hh

#include <protocols/splice/SpliceInTail.fwd.hh>
#include <protocols/splice/SpliceOut.hh>
#include <protocols/splice/SpliceOutTail.hh>
#include <protocols/splice/SpliceIn.hh>
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
#include <core/sequence/SequenceProfile.hh>
#include <protocols/task_operations/SeqprofConsensusOperation.fwd.hh>
#include <basic/database/open.hh>
#include <core/pack/task/PackerTask.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <protocols/moves/Mover.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>



namespace protocols {
namespace splice {



class SpliceInTail : virtual public protocols::splice::SpliceIn
{

public:
	SpliceInTail();
	void apply( Pose & pose ) override;
	virtual std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new SpliceInTail ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )override;
	virtual ~SpliceInTail();
	utility::vector1<core::Size>::const_iterator dbase_begin() const {return tail_dbase_subset_.begin();}
	utility::vector1<core::Size>::const_iterator dbase_end() const {return tail_dbase_subset_.end();}
	//void read_torsion_database();
	void set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc) override;
	void set_fold_tree_nodes(core::pose::Pose const & pose) override;
	core::Size set_anchor_res()override;
	//void minimize_segment(core::pose::Pose & pose);
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string mover_name();
	virtual core::Size find_dbase_entry(core::pose::Pose const & pose) override;
	virtual void assign_from_res_to_res(core::pose::Pose const & pose) override;
	void build_ideal_segment(core::pose::Pose & pose) override;



private:
	utility::vector1< core::Size > tail_dbase_subset_; //dflt false;
	DataccacheBoolDataOP end_tail_dbase_subset_;

};

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_Splice_hh
