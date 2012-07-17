// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh
/// @brief header file for CcdLoopClosureMover protocol
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_CcdLoopClosureMover_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_CcdLoopClosureMover_hh

// Unit Headers
#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.fwd.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/loops/Loops.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/moves/Mover.hh>

// ObjexxFCL Headers

// Utility headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
// AUTO-REMOVED #include <vector>

#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

/// Move these forward declarations to CcdLoopClosureMover.fwd.hh
class CcdLoopClosureMover;
typedef utility::pointer::owning_ptr< CcdLoopClosureMover > CcdLoopClosureMoverOP;
typedef utility::pointer::owning_ptr< CcdLoopClosureMover const > CcdLoopClosureMoverCOP;

//wrapper for fast_ccd_loop_closure
class CcdLoopClosureMover : public moves::Mover {
public:
  CcdLoopClosureMover(
     Loop const& loop_def,
     core::kinematics::MoveMapCOP mm
  );

	~CcdLoopClosureMover();

	void set_max_rama_score_increase( core::Real input_max_rama_score_increase );

	void set_max_total_delta(
	  core::Real input_max_delta_helix,
		core::Real input_max_delta_strand,
		core::Real input_max_delta_loop );

	void set_tolerance( core::Real input_tolerance );

	void set_ccd_cycles( core::Size input_ccd_cycles );

	void set_bRama_check( bool input_bRama_check );

	virtual void apply( core::pose::Pose &pose );
	virtual std::string get_name() const;

	core::Real forward_deviation() const {
		return forward_deviation_;
	}

	core::Real backward_deviation() const {
		return backward_deviation_;
	}

	bool success() const {
		return ( forward_deviation() < tolerance_ ) && ( backward_deviation() < tolerance_ );
	}

	Loop get_loop() const { return loop_ ;}
	core::Size get_loop_start() const { return loop_.start(); }
	core::Size get_loop_stop() const { return loop_.stop(); }
	core::Size get_loop_cut() const { return loop_.cut(); }
	core::Size get_loop_size() const { return loop_.size(); }
	core::Real get_loop_skip_rate() const { return loop_.skip_rate(); }
	std::string get_loop_extended() const { return ( (loop_.is_extended()) ? ("true") : ("false") );}
	core::Real get_max_rama_score_increase() const;
	core::Real get_tolerance() const;
	core::Real get_max_total_delta( std::string secstr ) const;
	core::Size get_ccd_cycles() const { return ccd_cycles_; }
	std::string bRama_check() const { return ( (bRama_check_) ? ("true") : ("false") ); }
	friend std::ostream &operator<< ( std::ostream &os, CcdLoopClosureMover const &mover );

private:
  Loop loop_;
  core::kinematics::MoveMapCOP movemap_;
  core::Real max_rama_score_increase_;
  core::Real max_total_delta_helix_;
  core::Real max_total_delta_strand_;
  core::Real max_total_delta_loop_;
  core::Real tolerance_;

  Size ccd_cycles_;
  bool bRama_check_;

  core::Real forward_deviation_; // output
  core::Real backward_deviation_; // output
  core::Real torsion_delta_;
  core::Real rama_delta_;

  Size actual_cycles_;
};

class CcdMover;
typedef utility::pointer::owning_ptr< CcdMover > CcdMoverOP;
typedef utility::pointer::owning_ptr< CcdMover const > CcdMoverCOP;

class CcdMover : public moves::Mover {
public:
	CcdMover(
		Loop const& loop_def,
		core::kinematics::MoveMapCOP mm
	);

	~CcdMover();

	virtual void apply( core::pose::Pose &pose );
	virtual std::string get_name() const;

private:
	core::Size total_moves_;
	Loop loop_;
	core::kinematics::MoveMapCOP movemap_;
};

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_CcdLoopClosureMover_hh
