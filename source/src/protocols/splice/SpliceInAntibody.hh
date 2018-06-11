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

#ifndef INCLUDED_protocols_splice_SpliceInAntibody_hh
#define INCLUDED_protocols_splice_SpliceInAntibody_hh

#include <protocols/splice/SpliceInAntibody.fwd.hh>
#include <protocols/splice/SpliceInTail.hh>
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



class SpliceInAntibody : virtual public protocols::splice::SpliceInTail
{

public:
	SpliceInAntibody();
	void apply( Pose & pose ) override;
	virtual std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new SpliceInAntibody ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	virtual ~SpliceInAntibody();
	core::Size find_dbase_entry(core::pose::Pose const & pose) override;
	core::Size database_entry()const {return database_entry_; }
	void database_entry( core::Size const d ){ database_entry_ = d; }
	//void load_from_checkpoint();
	utility::vector1<core::Size>::const_iterator dbase_begin() const {return dbase_subset_.begin();}
	utility::vector1<core::Size>::const_iterator dbase_end() const {return dbase_subset_.end();}
	bool dbase_iterate() const { return dbase_iterate_; }
	void dbase_iterate( bool const d ){ dbase_iterate_ = d; }
	bool min_seg() const { return min_seg_; }
	void min_seg( bool const d ){ min_seg_ = d; }
	//void read_torsion_database();
	void database_pdb_entry( std::string const & s ){ database_pdb_entry_ = s; }
	std::string database_pdb_entry() const { return database_pdb_entry_; }
	void scorefxn(core::scoring::ScoreFunctionOP sf) {scorefxn_ = sf;}
	core::scoring::ScoreFunctionOP scorefxn() const {return scorefxn_;}
	void set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc) override;
	void set_fold_tree_nodes(core::pose::Pose const & pose) override;
	core::Size set_anchor_res() override;
	void find_vl_vh_cut(core::pose::Pose pose);
	void update_vl_vh_cut();
	void find_disulfide_postions(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos);
	virtual void assign_from_res_to_res(core::pose::Pose const & pose) override;
	//void copy_stretch(core::pose::Pose & pose);
	void adjust_template_jump(core::pose::Pose & pose);
	void build_ideal_segment(core::pose::Pose & pose) override;
	void antibody_DB(std::string const & s){antibody_DB_=s;}
	std::string antibody_DB(){return antibody_DB_;}
	void adjust_n_ter_tail_length(core::pose::Pose & pose);


	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string mover_name();


private:
	core::Size vl_vh_cut_=0;
	utility::vector1<core::Size> pose_cys_pos_,template_cys_pos_;
	std::string antibody_DB_;
	int tail_diff_=0; //if source antibody has a tail that is shorter than that of 2BRR (the template AB) we need to acount for that)


};

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_Splice_hh
