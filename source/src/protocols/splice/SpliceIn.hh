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
/// @brief Takes a conformtion phi/spi angles from an input db file and applies to the corresponding segment in the pose

#ifndef INCLUDED_protocols_splice_SpliceIn_hh
#define INCLUDED_protocols_splice_SpliceIn_hh

#include <protocols/splice/SpliceIn.fwd.hh>
#include <protocols/splice/SpliceOut.hh>
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
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <protocols/moves/Mover.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>



namespace protocols {
namespace splice {



class SpliceIn :  virtual public protocols::splice::SpliceOut
{

public:
	SpliceIn();
	void apply( Pose & pose ) override;
	virtual std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new SpliceIn ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	virtual ~SpliceIn();
	virtual core::Size find_dbase_entry(core::pose::Pose const & pose); //get pdb entry from database
	virtual void assign_from_res_to_res(core::pose::Pose const & pose); //assign the start and end residues of the Spliced in segment
	core::Size database_entry()const {return database_entry_; } //setter for db entry
	void database_entry( core::Size const d ){ database_entry_ = d; } //getter for db entry
	void load_from_checkpoint(); //restart run from checkpoint
	utility::vector1<core::Size>::const_iterator dbase_begin() const {return dbase_subset_.begin();} //
	utility::vector1<core::Size>::const_iterator dbase_end() const {return dbase_subset_.end();}
	bool dbase_iterate() const { return dbase_iterate_; }
	void dbase_iterate( bool const d ){ dbase_iterate_ = d; }
	bool min_seg() const { return min_seg_; }
	void min_seg( bool const d ){ min_seg_ = d; }
	virtual void minimize_segment(core::pose::Pose & pose); //minimize segment after splice in
	virtual void read_torsion_database();
	void database_pdb_entry( std::string const & s ){ database_pdb_entry_ = s; }
	std::string database_pdb_entry() const { return database_pdb_entry_; }
	void scorefxn(core::scoring::ScoreFunctionOP sf) {scorefxn_ = sf;}
	core::scoring::ScoreFunctionOP scorefxn() const {return scorefxn_;}
	void set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc) override; //apply loop length change before changing torsion angles
	void set_fold_tree_nodes(core::pose::Pose const & pose) override; //how to build the fold tree
	virtual core::Size set_anchor_res() override; //set anchor res for coordinate constraints

	virtual void build_ideal_segment(core::pose::Pose & pose) override; //build segment using ideal bond length and angles
	void rtmin( core::pose::Pose & pose,core::pack::task::TaskFactoryOP tf); //apply rotamer trial min mover
	void pack( core::pose::Pose & pose,core::pack::task::TaskFactoryOP tf); //apply packing


	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string mover_name();


protected:
	utility::vector1< core::Size > dbase_subset_;
	utility::vector1< core::Size >::const_iterator current_dbase_entry_;
	bool dbase_iterate_; //dflt false;
	DataccacheBoolDataOP end_dbase_subset_;
	std::string checkpointing_file_;
	bool min_seg_;// dflt true, if set to true
	utility::vector1< ResidueBBDofs > torsion_database_;
	bool first_pass_; // dflt true;
	core::Size database_entry_; //dflt 0; in which case tests a random entry in each apply
	std::string database_pdb_entry_; // dflt ""; e.g., "1yihl" specify this only if you want just one loop to be spliced
	utility::vector1< ResidueBBDofs > tail_torsion_database_;
	utility::vector1< int > delta_lengths_; // dflt empty; change loop length by how much? 0 is always assumed
	core::Size allowed_cuts_; //dflt=1. if the number of chain_breaks in the psoe is more than the allowed number then fail

};

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_Splice_hh
