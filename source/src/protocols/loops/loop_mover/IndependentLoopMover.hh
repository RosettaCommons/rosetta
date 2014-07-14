// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/IndependentLoopMover.hh
/// @brief
/// @author Mike Tyka
/// @author James Thompson

#ifndef INCLUDED_protocols_loops_loop_mover_IndependentLoopMover_hh
#define INCLUDED_protocols_loops_loop_mover_IndependentLoopMover_hh

#include <protocols/loops/loop_mover/IndependentLoopMover.fwd.hh>

#include <protocols/loops/loop_mover/LoopMover.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/LoopsFileIO.fwd.hh>

#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {
namespace loop_mover {

// A subclass where all the loops are modelled independently from each
// other.  Deriving classes should not overload apply but instead
// model_loop. This is an unfortunate design decision, because there's no
// way to guarantee that deriving classes do the right thing and don't
// override the apply method.
class IndependentLoopMover : public LoopMover {
public:
	IndependentLoopMover();
	IndependentLoopMover( LoopsOP loops_in );
  IndependentLoopMover( utility::vector1< bool > const& selection );
	IndependentLoopMover( LoopsFileData const & lfd );
	IndependentLoopMover( GuardedLoopsFromFileOP guarded_loops );
	
	//destructor
	virtual ~IndependentLoopMover();
	
	void set_defaults();

	/// @brief Apply the loop-build protocol to the input pose
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	/// Accessors:

	void set_build_attempts_( int value ) { build_attempts_ =  value; }
	void set_grow_attempts_( int value )  { grow_attempts_ =  value; }
	void set_accept_aborted_loops_( bool value ) {   accept_aborted_loops_ =  value; }
	void set_strict_loops( bool value ) { strict_loops_ =  value; }
	void set_random_order_( bool value ) { random_order_ =  value; }
	void set_build_all_loops_( bool value ) {   build_all_loops_ =  value; }
	void set_loop_combine_rate_( core::Real value ) {   loop_combine_rate_ =  value; }

	int  get_build_attempts() const { return build_attempts_; }
	int  get_grow_attempts() const { return grow_attempts_; }
	bool get_accept_aborted_loops() const { return accept_aborted_loops_; }
	bool get_strict_loops() const { return strict_loops_; }
	bool get_random_order() const { return random_order_; }
	bool get_build_all_loops() const { return build_all_loops_; }
	bool get_loop_combine_rate() const { return loop_combine_rate_; }
	bool get_all_loops_closed() const {return all_loops_closed_; }

private:

	/// select loops to be built
	void select_loops( Loops & selected_loops );

	/// Try loopbuilding n times before extending
	int build_attempts_;

	/// Try extending n times
	int grow_attempts_;

	/// danger - can lead to infinite loops !
	bool accept_aborted_loops_;

	/// Grow loops outwards if building fails.
	bool strict_loops_;

	/// Randomise loop build order
	bool random_order_;

	/// Force to build all loops (i.e. ignore skiprate)
	bool build_all_loops_;

	/// determines if all the loops were closed in the last time mover was
	/// called.
	bool all_loops_closed_;
	/// Loop combine rate
	core::Real loop_combine_rate_;
    
protected:
	virtual
	LoopResult
	model_loop(
		core::pose::Pose & pose,
		protocols::loops::Loop const & loop
	) = 0;

};

} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_IndependentLoopMover_hh
