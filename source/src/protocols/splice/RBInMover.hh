// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RBInMover.hh
/// @brief Get rigid body orientations vH antibody chains when vL is aligned.

#ifndef INCLUDED_protocols_splice_RBInMover_hh
#define INCLUDED_protocols_splice_RBInMover_hh

#include <protocols/splice/RBInMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <algorithm>
#include <core/kinematics/Jump.hh>

// C++ Headers
namespace protocols {
namespace splice {

class RBInMover : public protocols::moves::Mover {
public:
	RBInMover();
	~RBInMover() override;

	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string RBDatabase_fname() const {return RB_dbase_; }
	void RBDatabase_fname( std::string const & s ){ RB_dbase_ = s; }

	void RB_dbase( std::string const & s ){ RB_dbase_ = s; }
	std::string RB_dbase() const{ return RB_dbase_; }
	core::Size from_entry() const{ return from_entry_; }
	void from_entry( core::Size const s ){ from_entry_ = s; }
	void to_entry( core::Size const s ){ to_entry_ = s ;}
	core::Size to_entry() const{ return to_entry_;}
	bool randomize() const {return randomize_; }
	void randomize( bool const b ){ randomize_ = b;}
	utility::vector1< std::pair <std::string,core::kinematics::Jump> > jump_library() const;
	void jump_library( utility::vector1< std::pair <std::string,core::kinematics::Jump> > j );
	void checkpointing_file( std::string const & s ){ checkpointing_file_ = s; }
	std::string checkpointing_file() const{ return checkpointing_file_; }

	void set_fold_tree( core::pose::Pose & pose ) ;

	void modify_foldtree( bool const m ){ modify_foldtree_ = m; }
	bool modify_foldtree() const { return modify_foldtree_; }
	void find_disulfide_postions(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos);
	core::Size find_vl_vh_cut(core::pose::Pose pose) const;
	utility::vector1<std::array<int, 3>> set_fold_tree_nodes(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos, core::Size vl_vh_cut);

	void db_entry( std::string s) {db_entry_=s ;}

	std::string db_entry(){ return db_entry_; }
	void select_entry_by_name();


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void init(); /// sets the entry order; if the jump_library is already populated returns without doing anything
	bool checkpoint_recovery(); /// recover from checkpointing. If checkpointing is off, does nothing. Returns true if recovered from checkpoint
	bool checkpoint() const; // dump the checkpoint information to disk. Return whether or not checkpointing took place

	std::string RB_dbase_; //dflt ""
	core::Size from_entry_, to_entry_, current_entry_; //dflt 1,1,1
	bool randomize_; //dflt true
	utility::vector1< std::pair <std::string,core::kinematics::Jump> > jump_library_; //dflt empty;


	std::string checkpointing_file_; //dflt "" in which case we're not checkpointing
	bool modify_foldtree_; //dflt true; use the existing fold tree or modify it to go by the variable domain disulfides.
	std::string db_entry_;
};


} // simple_moves
} // protocols

#endif
