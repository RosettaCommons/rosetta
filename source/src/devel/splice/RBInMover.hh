// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RBInMover.hh
/// @brief Get rigid body orientations vH antibody chains when vL is aligned.

#ifndef INCLUDED_devel_splice_RBInMover_hh
#define INCLUDED_devel_splice_RBInMover_hh

#include <devel/splice/RBInMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <algorithm>
#include <core/kinematics/Jump.hh>

// C++ Headers
namespace devel {
namespace splice {
    
class RBInMover : public protocols::moves::Mover {
public:
	RBInMover();
    ~RBInMover();

    
    //	RBInMover(core::Real const min_in , core::Real const max_in);
    
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

    std::string RBDatabase_fname() const {return RB_dbase_; }
    void RBDatabase_fname( std::string const s ){ RB_dbase_ = s; }
   
		void RB_dbase( std::string const s ){ RB_dbase_ = s; }
		std::string RB_dbase() const{ return RB_dbase_; }
		core::Size from_entry() const{ return from_entry_; }
		void from_entry( core::Size const s ){ from_entry_ = s; }
		void to_entry( core::Size const s ){ to_entry_ = s ;}
 		core::Size to_entry() const{ return to_entry_;}
 		bool randomize() const {return randomize_; }
		void randomize( bool const b ){ randomize_ = b;}
    //void jump_dbase_fname( std::string const s ) { jump_dbase_fname_ = s ; }
    //std::string jump_dbase_fname() const{ return jump_dbase_fname_; }
	  void init(); /// sets the entry order
	  utility::vector1< core::kinematics::Jump > jump_library() const;
    void jump_library( utility::vector1< core::kinematics::Jump > j );
private:
    std::string RB_dbase_; //dflt ""
		core::Size from_entry_, to_entry_, current_entry_; //dflt 1,1,1
		bool randomize_; //dflt true
    //std::string jump_dbase_fname_;
	  utility::vector1< core::kinematics::Jump > jump_library_; //dflt empty;
};


} // simple_moves
} // protocols

#endif
