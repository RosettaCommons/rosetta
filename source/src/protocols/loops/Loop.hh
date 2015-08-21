// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/Loop.hh
/// @brief
/// @author Chu Wang
/// @author Mike Tyka
/// @author James Thompson

#ifndef INCLUDED_protocols_loops_Loop_HH
#define INCLUDED_protocols_loops_Loop_HH

// Unit headers
#include <protocols/loops/Loop.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#ifdef WIN32
#include <functional>
#endif

#include <iostream>
//#include <ostream>

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {

/// @brief Bare-bones representation of a loop
struct SerializedLoop {
	core::Size start;
	core::Size stop;
	core::Size cut;
	core::Real skip_rate;
	bool extended;
};

/// single loop definition
class Loop : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Loop();
	/// default constructor
	Loop():
		start_(0),
		stop_(0),
		cut_(0),
		skip_rate_( 0.0 ),
		extended_(false)
	{}

	Loop( SerializedLoop loop ):
		start_( loop.start ),
		stop_( loop.stop ),
		cut_( loop.cut ),
		skip_rate_( loop.skip_rate ),
		extended_( loop.extended )
	{}

	/// input constructor
	Loop(
		core::Size const start_in, core::Size const stop_in,
		core::Size const cut_in = 0, core::Real skip_rate = 0.0,
		bool const extended_in = false
	):
		start_( start_in ),
		stop_( stop_in ),
		cut_( cut_in ),
		skip_rate_( skip_rate ),
		extended_( extended_in )
	{}

	inline bool is_extended() const { return extended_; }
	inline core::Size start() const { return start_; }
	inline core::Size stop() const { return stop_; }
	inline core::Size cut() const { return cut_; }
	inline core::Size size() const { return stop_ - start_ + 1; }
	inline core::Real skip_rate() const { return skip_rate_; }

	inline void set_extended( bool input    ) { extended_ = input; }
	inline void set_start( core::Size input ) { start_ = input; }
	inline void set_stop( core::Size input  ) { stop_ = input; }
	inline void set_cut( core::Size input   ) { cut_ = input; }

	/// @brief Assuming that the loop represents a contiguous stretch of residues,
	/// returns the length. Makes no assumptions about directionality. That is,
	/// Loop(3,8).length() == Loop(8,3).length(). Constant time procedure.
	core::Size length() const {
		Size m = std::min(start(), stop());
		Size n = std::max(start(), stop());
		return n - m + 1;
	}

	/// @brief Returns true if the loop's elements are increasing
	bool increasing() const {
		return start() <= stop();
	}

	/// @brief Returns true if the loop's elements are decreasing
	bool decreasing() const {
		return !increasing();
	}

	/// @brief Returns the midpoint of the loop
	core::Size midpoint() const {
		return increasing() ? start() + length() / 2 : stop() - length() / 2;
	}

	bool operator< ( Loop const& larger ) const {
		return ( size() < larger.size() ? true :
			( start() < larger.start() ? true : ( cut() < larger.cut() ) ) );
	}
	bool operator==( Loop const& other ) const {
		return ( size() == other.size() &&
			start() == other.start() && cut() == other.cut() );
	}
	bool operator!=( Loop const& other ) const {
		return !(*this == other );
	}

	/// @brief add all residues within this loop definition into selection
	void get_residues( utility::vector1< Size>& selection ) const;

	/// @brief switch DOF_Type for residues in loop. id::CHI, id::BB --- don't
	/// use with id::JUMP
	void switch_movemap( core::kinematics::MoveMap& movemap, core::id::TorsionType, bool allow_moves = true ) const;

	// @brief Autochoose a cutpoint using the secondary structure of the pose.
	void choose_cutpoint( core::pose::Pose const & pose );

	/// @brief Autochoose a cutpoint using the secondary structure of the pose
	/// unless cutpoint is already set
	void auto_choose_cutpoint( core::pose::Pose const & pose ){
		if ( cut_ <= 0 ) choose_cutpoint ( pose );
	}

	bool is_terminal( core::pose::Pose const & pose ) const;

	/// @brief  Generate a string representation of Loop
	virtual void show( std::ostream & output=std::cout ) const;

private:
	core::Size start_;
	core::Size stop_;
	core::Size cut_;
	core::Real skip_rate_;
	bool extended_;
}; // Loop

//////////////////////////////////////////////////////////////////////
std::ostream & operator<<( std::ostream & os, Loop const & loop );

/// @brief Orders loops by start position
class RationalLoopComparator : public std::binary_function<double, double, bool> {
public:
	bool operator()(Loop x, Loop y) {
		return x.start() < y.start();
	}
};


/// @brief used to sort Loops by start-res
class Loop_lt : public std::binary_function<double, double, bool> {
public:
	bool operator()(Loop x, Loop y) {
		return (x.start() < y.stop());  // so wrong...
	}
};

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_Loop_HH
