// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RBOutMover.hh
/// @brief Get rigid body orientations vH antibody chains when vL is aligned.

#ifndef INCLUDED_devel_splice_RBOutMover_hh
#define INCLUDED_devel_splice_RBOutMover_hh

#include <devel/splice/RBOutMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <algorithm>

// C++ Headers
namespace devel {
namespace splice {
    
class RBOutMover : public protocols::moves::Mover {
public:
	RBOutMover();
    ~RBOutMover();

    
    //	RBOutMover(core::Real const min_in , core::Real const max_in);

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
    
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

    std::string template_pdb_fname() const {return template_pdb_fname_; }
    void template_pdb_fname( std::string const s ){ template_pdb_fname_ = s; }
    
    void jump_dbase_fname( std::string const s ) { jump_dbase_fname_ = s ; }
    std::string jump_dbase_fname() const{ return jump_dbase_fname_; }
private:
    std::string template_pdb_fname_;
    std::string jump_dbase_fname_;
    //	core::Real min_value_,max_value_; // dflt -1, 5
};

/// utility function for finding disulfide bonded pairs
utility::vector1< std::pair< core::Size, core::Size > >
find_disulfs_in_range( core::pose::Pose const & pose, core::Size const start, core::Size const end );

} // simple_moves
} // protocols

#endif
