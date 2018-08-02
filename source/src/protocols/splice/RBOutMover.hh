// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RBOutMover.hh
/// @brief Get rigid body orientations vH antibody chains when vL is aligned.

#ifndef INCLUDED_protocols_splice_RBOutMover_hh
#define INCLUDED_protocols_splice_RBOutMover_hh

#include <protocols/splice/RBOutMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <algorithm>
#include <core/kinematics/Jump.hh>
#include <basic/datacache/DataMap.hh>

// C++ Headers
namespace protocols {
namespace splice {

class RBOutMover : public protocols::moves::Mover {
public:
	RBOutMover();
	~RBOutMover() override;

	core::kinematics::Jump get_disulf_jump( Pose & pose, core::pose::Pose & template_pose );


	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string template_pdb_fname() const {return template_pdb_fname_; }
	void template_pdb_fname( std::string const & s ){ template_pdb_fname_ = s; }

	void jump_dbase_fname( std::string const & s ) { jump_dbase_fname_ = s ; }
	std::string jump_dbase_fname() const{ return jump_dbase_fname_; }

	bool jump_from_foldtree() const{ return jump_from_foldtree_; }
	void jump_from_foldtree( bool const jfft ){ jump_from_foldtree_ = jfft;}
	void find_disulfide_postions(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos);
	core::Size find_vl_vh_cut(core::pose::Pose pose);
	utility::vector1<std::array<int, 3>> set_fold_tree_nodes(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos, core::Size vl_vh_cut);
	void superimpose_source_on_pose( core::pose::Pose const & target_pose,core::Size target_from_res,core::Size target_to_res, core::pose::Pose & template_pose ,
		core::Size template_from_res,core::Size template_to_res);
	std::string
	get_name() const override;
	bool debug(){return debug_;}
	void debug(bool b){debug_=b;}

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string template_pdb_fname_ = "";
	std::string jump_dbase_fname_ = "";
	bool jump_from_foldtree_ = false; /// if true, extract the jump defined by the fold tree rather than imposing a new jump
	// core::Real min_value_,max_value_; // dflt -1, 5
	bool debug_ = false;
};

/// utility function for finding disulfide bonded pairs
utility::vector1< std::pair< core::Size, core::Size > >
find_disulfs_in_range( core::pose::Pose const & pose, core::Size const start, core::Size const end );

} // splice
} // protocols

#endif
