// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_BoolMover_hh
#define INCLUDED_protocols_moves_BoolMover_hh

// Unit Headers
// AUTO-REMOVED #include <protocols/moves/Mover.fwd.hh>

// Package headers
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/MoverStatistics.hh>

// Project headers
// AUTO-REMOVED #include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>

// ObjexxFCL Headers

// C++ Headers
// AUTO-REMOVED #include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

//Auto Headers
#include <iostream>


namespace protocols {
namespace moves {

//temporary code duplication of Mover class to return bool from apply..
// want to make this interface the one and only Mover at some point
//  This is completely obsolet now... the new jd2 design takes care of this issue.
// and exceptions of course
class _BoolMover : public utility::pointer::ReferenceCount {
public:
	_BoolMover();
	virtual ~_BoolMover();

	_BoolMover( std::string const & type );

	_BoolMover( _BoolMover const & other );

	virtual bool apply( core::pose::Pose & ) = 0;

	std::string const & type() const { return type_; }

	/// @brief A tag is a unique identifier used to identify structures produced
	/// by this Mover. get_current_tag() returns the tag, and set_current_tag( std::string tag )
	/// sets the tag.
	std::string const & get_current_tag() const { return current_tag_; }
	void set_current_tag( const std::string& new_tag ) { current_tag_ = new_tag; }

	///@brief setters and getters for poses contained for rms
	void set_input_pose( core::pose::PoseCOP pose );
	void set_native_pose( core::pose::PoseCOP pose );



	core::pose::PoseCOP get_input_pose() const;
	core::pose::PoseCOP get_native_pose() const;

	///@brief: Unit test support function.  Apply one move to a given pose.
	///  			 Allows extra test specific functions to be called before applying
  virtual void test_move( core::pose::Pose & pose ) {
		apply( pose );
	}

	//protected:
	//OL 4/23/08 made this public. it is not really a safety issue to have that
	//protected and it allows more detail in MC diagnosis
	void type( const std::string & type_in ) { type_ = type_in; }

private:
	std::string type_; // should be const
	std::string current_tag_;

	core::pose::PoseCOP input_pose_;
	core::pose::PoseCOP native_pose_;
}; // end Mover base class

}
}

#endif
