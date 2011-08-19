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
/// @author Mike Tyka

#ifndef INCLUDED_protocols_loops_LoopMover_hh
#define INCLUDED_protocols_loops_LoopMover_hh

#include <protocols/loops/LoopMover.fwd.hh>

#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>  // Changing from .fwd.hh --> .hh to make python bindings compile
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <protocols/checkpoint/CheckPointer.hh>

// C++ Headers


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {


enum LoopResult { Success, CriticalFailure, Failure, ExtendFailure };


/// @brief The loop-rebuild protocol
class LoopMover: public protocols::moves::Mover {
public: // typedefs

	typedef core::kinematics::MoveMap MoveMap;

public:

	LoopMover(
		protocols::loops::Loops loops_in
	) : Mover(),
		loops_( loops_in ),
		checkpoints_("LoopMover"),
		loops_from_observer_cache_(false)
	{
		Mover::type("LoopMover");
	}

	void set_scorefxn( const core::scoring::ScoreFunctionOP score_in ) {
		scorefxn_ = score_in;
	}

	/// @brief Apply the loop-build protocol to the input pose
	void apply( core::pose::Pose & ) {}

	virtual std::string get_name() const;

	const protocols::loops::Loops & loops() const {
		return loops_;
	}

	/// @brief Extend a loop
	virtual void set_extended_torsions(
		core::pose::Pose & pose,
		Loop const & loop
	);

public: // fragment libraries

	/// @brief add a fragment set
	void add_fragments( core::fragment::FragSetOP fragset );

	/// @brief clear the list of fragment sets
	void clear_fragments();

public: // movemap management

	/// @brief <b>explicit</b> False settings in this MoveMap will override any
	///  automatically generated MoveMap settings during the loop modeling
	///  protocol
	inline
	MoveMap const & false_movemap() const {
		return false_movemap_;
	}

	/// @brief <b>explicit</b> False settings in this MoveMap will override any
	///  automatically generated MoveMap settings during the loop modeling
	///  protocol
	inline
	void false_movemap( MoveMap const & mm ) {
		false_movemap_ = mm;
	}

public: // checkpointing

	checkpoint::CheckPointer & get_checkpoints() {
		return checkpoints_;
	}

protected: // movemap management

	/// @brief import the false_movemap's <b>explicit</b> False settings into the
	///  given MoveMap
	/// @return The number of False settings imported.
	inline
	Size enforce_false_movemap( MoveMap & mm ) const {
		return mm.import_false( false_movemap_ );
	}

	/// @author flo, march 2011
	/// @brief allow the loops to be set from the segments
	/// stored in the poses observer cache. makes it possible
	/// to have LoopMovers be part of parser protocols
	/// where the loops were determined by some previous on the
	/// fly step
	void
	set_loops_from_pose_observer_cache( core::pose::Pose const & pose );

protected: // data

	protocols::loops::Loops loops_;
	core::scoring::ScoreFunctionOP scorefxn_;

	checkpoint::CheckPointer checkpoints_;

	// it's really dumb that this is a std::vector rather a vector1
	std::vector< core::fragment::FragSetOP > frag_libs_;

	bool loops_from_observer_cache_;

private: // data

	/// @brief <b>explicit</b> False settings in this MoveMap will override any
	///  automatically generated MoveMap settings during the loop modeling
	///  protocol
	MoveMap false_movemap_;

}; // class LoopMover

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_LoopMover_HH
