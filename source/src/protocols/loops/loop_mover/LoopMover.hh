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

#ifndef INCLUDED_protocols_loops_loop_mover_LoopMover_hh
#define INCLUDED_protocols_loops_loop_mover_LoopMover_hh

#include <protocols/loops/loop_mover/LoopMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/conformation/ppo_torsion_bin.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/LoopsFileIO.fwd.hh>

#include <protocols/checkpoint/CheckPointer.fwd.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

#ifdef WIN32
//#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#endif

// Utility headers
#include <utility/vector1.hh>

// C++ Headers


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {
namespace loop_mover {

enum LoopResult { Success, CriticalFailure, Failure, ExtendFailure };


/// @brief The loop-rebuild protocol
class LoopMover: public protocols::moves::Mover {
public: // typedefs

	typedef core::kinematics::MoveMapOP MoveMapOP;

public:

	LoopMover();
	LoopMover( protocols::loops::LoopsOP loops_in );
	LoopMover( protocols::loops::LoopsFileData const & loops_from_file );
	LoopMover( protocols::loops::GuardedLoopsFromFileOP guarded_loops );

	/// @brief copy ctor
	LoopMover( LoopMover const & rhs );

	/// @brief assignment operator
	LoopMover & operator=( LoopMover const & rhs );

	//destructor
	virtual ~LoopMover();

	virtual std::string get_name() const;

	/// @brief Inform the GuardedLoopsFromFile object that it is not in charge of
	/// updating its Loops object at the beginning of apply()
	void set_guarded_loops_not_in_charge();

	/// @brief Apply the loop-build protocol to the input pose
	virtual void apply( core::pose::Pose & ) = 0;

	void set_scorefxn( const core::scoring::ScoreFunctionOP score_in );
	const core::scoring::ScoreFunctionOP & scorefxn() const;

	/// @brief Set the loops pointer by giving the LoopMover resolved loop indices; implicity sets
	/// the GuardedLoopsFromFile object into a "not in charge" state (since something else
	/// must be controlling the the Loops object).  The GuardedLoopFromFile object copies the pointer,
	/// not the data.
	void loops( protocols::loops::LoopsOP lptr );

	/// @brief Set the loops by giving the LoopMover unresolved loop indices (which cannot be resolved until apply() ).
	void loops( LoopsFileData const & loop_file_data );

	/// @brief Set the guarded_loops pointer
	void loops( protocols::loops::GuardedLoopsFromFileOP guarded_loops );

	/// @brief Accessor for the loops data.  Requires that the loop indices have been resolved; do not call this before
	/// apply() has been called.
	protocols::loops::LoopsCOP loops() const;

	/// @brief non-const accessor for the loops data.  Requires that the loop indices have been resolved; do not call this
	/// before apply() has been called.
	protocols::loops::LoopsOP loops();

	const utility::vector1< core::fragment::FragSetOP > & frag_libs() const;

	/// @brief create a string representing the torsion bins (ABEGX) for the loops defined in the
	/// guarded_loops_ object.  A sentinal value of ppo_torbin_U is used to mark the boundaries
	/// between the separate loops (unused as long as there's only one Loop object in the
	/// guarded_loops_ object).
	core::conformation::torsion_bin_string
	torsion_features_string( core::pose::Pose const & pose ) const; // AS

	/// @brief Extend a loop
	virtual void set_extended_torsions( core::pose::Pose & pose, Loop const & loop );

public: // fragment libraries

	/// @brief add a fragment set
	void add_fragments( core::fragment::FragSetOP fragset );

	/// @brief clear the list of fragment sets
	void clear_fragments();

public: // movemap management

	/// @brief <b>explicit</b> False settings in this MoveMap will override any
	///  automatically generated MoveMap settings during the loop modeling
	///  protocol
	MoveMapOP const & false_movemap() const;

	/// @brief <b>explicit</b> False settings in this MoveMap will override any
	///  automatically generated MoveMap settings during the loop modeling
	///  protocol
	void false_movemap( MoveMapOP const & mm );

public: // checkpointing
	checkpoint::CheckPointerOP & get_checkpoints();

protected: // movemap management

	/// @brief import the false_movemap's <b>explicit</b> False settings into the
	///  given MoveMap
	/// @return The number of False settings imported.
	Size enforce_false_movemap( MoveMapOP & mm ) const;

	/// @author flo, march 2011
	/// @brief allow the loops to be set from the segments
	/// stored in the poses observer cache. makes it possible
	/// to have LoopMovers be part of parser protocols
	/// where the loops were determined by some previous on the
	/// fly step
	void set_loops_from_pose_observer_cache( core::pose::Pose const & pose );

	bool use_loops_from_observer_cache() const;
	void set_use_loops_from_observer_cache( bool const loops_from_observer_cache );

	virtual basic::Tracer & tr() const = 0;


	/// @brief Turn the unresolved loop indices read in from disk into pose-specific
	/// loop indices.  Must be called by derived classes at the beginning of apply.
	void resolve_loop_indices( core::pose::Pose const & p );

private:

	void init();
	void initForEqualOperatorAndCopyConstructor(LoopMover & lhs, LoopMover const & rhs);


private: // data

	GuardedLoopsFromFileOP guarded_loops_;

	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< core::fragment::FragSetOP > frag_libs_;
	checkpoint::CheckPointerOP checkpoints_;
	bool loops_from_observer_cache_;

	/// @brief <b>explicit</b> False settings in this MoveMap will override any
	///  automatically generated MoveMap settings during the loop modeling
	///  protocol
	MoveMapOP false_movemap_;

}; // class LoopMover

void
loops_set_chainbreak_weight( core::scoring::ScoreFunctionOP scorefxn, core::Size const round = 1 );

} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_LoopMover_HH
