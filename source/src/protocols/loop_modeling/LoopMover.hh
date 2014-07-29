// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopMover_HH
#define INCLUDED_protocols_loop_modeling_LoopMover_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <boost/utility.hpp>

// C++ headers
#include <string>

namespace protocols {
namespace loop_modeling {

/// @brief Represent the types of fold tree supported by LoopMover.
/// @details This enum supports the bitwise `and` and `or` operators.  Strictly 
/// speaking, this is an abuse of enum because both operators can return values 
/// that are not contained in the enum.  But this design meets two key 
/// requirements.  First, programmers have to work with the enumerated values.  
/// Second, bitwise logic can be used to easily combine flags and to determine 
/// how composite requests should be satisfied.

enum FoldTreeRequest {
	FTR_LOOPS_WITH_CUTS = 0x01,
	FTR_SIMPLE_TREE     = 0x02,
	FTR_DONT_CARE       = 0xFF
};

/// @brief Implement the bitwise `or` operator for fold tree requests.
inline FoldTreeRequest operator | (FoldTreeRequest a, FoldTreeRequest b) {
	return static_cast<FoldTreeRequest>(
			static_cast<Size>(a) | static_cast<Size>(b));
}

/// @brief Implement the bitwise `and` operator for fold tree requests.
inline FoldTreeRequest operator & (FoldTreeRequest a, FoldTreeRequest b) {
	return static_cast<FoldTreeRequest>(
			static_cast<Size>(a) & static_cast<Size>(b));
}

/// @brief Base class for loop-sampling algorithms.
/// @details Classes that inherit from LoopMover can plug into the LoopProtocol 
/// framework.  The only method that subclasses need to provide is do_apply().  
/// Note that the apply() method itself cannot be overwritten, because 
/// LoopMover uses it to do some useful setup and teardown work.  Instead, 
/// there are two versions of do_apply() that subclasses may implement.  Both 
/// return a boolean to indicate whether or not the move succeeded.  The first 
/// version of do_apply() accepts only a Pose and is expected to operate on all 
/// of the loops returned by get_loops().  The second version accepts a Pose 
/// and a Loop, and is only expected to operate on the given loop.  By default, 
/// the first version simply iterates through the loops provided by get_loops() 
/// and calls the second version on each one.  This means that if the first 
/// version is overwritten, the second version will no longer be called.  If 
/// neither method is reimplemented, a runtime error will be thrown.
///
/// LoopMover provides a handful of features that could be useful to a loop 
/// sampling algorithm.  As mentioned above, the get_loops() method returns the 
/// loops that should be sampled.  There are also a number of methods provided 
/// for controlling how the fold tree is setup up.  The request_fold_tree() 
/// method can be reimplemented to return an enum telling what kind of fold 
/// tree this mover requires.  When apply() is called for the first time, a 
/// compatible fold tree will be configured.  This behavior is disabled if 
/// trust_fold_tree() is called beforehand, in which case responsibility for 
/// constructing a compatible fold tree is passed to the calling code.

class LoopMover :
	public protocols::moves::Mover, protected boost::noncopyable {

public:

	/// @brief Default constructor.
	LoopMover();
	
	/// @brief Default destructor.
	~LoopMover();
	
	/// @brief Sample the pose in the regions specified by get_loops().
	/// @details The parent class apply() method automatically sets up a fold 
	/// tree (if necessary) and keeps track of whether or not the move succeeded.  
	/// Child classes should reimplement do_apply() instead of this method.
	void apply(Pose & pose);

	/// @brief Return the name of this mover.
	virtual string get_name() const { return "LoopMover"; }

public:

	/// @brief Return true if the previous move was successful.
	bool was_successful() const;
	
	/// @brief Return the loops to be sampled on the next call to apply().
	Loops get_loops() const;

	/// @brief Return the score function to be used on the next call to apply().
	core::scoring::ScoreFunctionCOP get_score_function() const;

	/// @brief Return the score function to be used on the next call to apply() 
	/// (non-const access).
	core::scoring::ScoreFunctionOP get_score_function();

	/// @brief Set the loop to be sampled on the next call to apply().
	void set_loop(Loop const & loop);

	/// @brief Set the loops to be sampled on the next call to apply().
	void set_loops(Loops const & loops);

	/// @brief Set the score function to be used on the next call to apply().
	void set_score_function(core::scoring::ScoreFunctionOP sfxn);

public:

	/// @brief Return an enum representing the kind of fold tree that is 
	/// compatible with this mover.
	/// @details The FoldTreeRequest enum values can be combined using the 
	/// bitwise logical operators.  For example, you can request either the 
	/// standard fold tree or a simple fold tree with `FTR_LOOPS_WITH_CUTS | 
	/// FTR_SIMPLE_TREE.`
	virtual FoldTreeRequest request_fold_tree() const;

	/// @brief Promise that the calling code will setup a fold tree compatible 
	/// with request_fold_tree().  If this method is not called, most movers will 
	/// setup a fold tree on their own the first time apply() is called.
	void trust_fold_tree();

	/// @brief Setup the given pose with a fold tree that is compatible with the 
	/// given loops and requests.
	static void setup_fold_tree(
			Pose & pose, Loops const & loops, FoldTreeRequest request);

protected:

	/// @brief Perform the loop sampling move.  This method can be overwritten in 
	/// child classes.
	virtual bool do_apply(Pose & pose);

	/// @brief Perform the loop sampling move.  This method can be overwritten in 
	/// child classes.
	virtual bool do_apply(Pose & pose, Loop const & loop);

	/// @brief Disable all the LoopMover score function related methods.
	/// @details This is useful in LoopMover subclasses like LoopModeler that 
	/// don't fit the "one score function" paradigm.  When score function 
	/// management is disabled, the get_score_function() and set_score_function() 
	/// methods print a warning to the tracer and not do anything else.  Nested 
	/// loop movers will need to be manually given score functions.
	void dont_manage_score_function();

	/// @brief Indicate that the given loop mover is used within do_apply().
	/// @details Once registered, the nested mover will be always be configured 
	/// with the same loops, the same score function, and the same fold tree as 
	/// this loop mover.  Fold tree requests made by nested movers will be taken 
	/// into account by the parent.  This synchronization is entirely managed by 
	/// LoopMover, so subclasses don't need to do anything special.
	template <class LoopMoverSubclassOP>
	LoopMoverSubclassOP register_nested_loop_mover(LoopMoverSubclassOP mover) {
		mover->trust_fold_tree();
		mover->set_loops(loops_);
		mover->set_score_function(score_function_);
		nested_movers_.push_back(mover);
		return mover;
	}

	/// @brief Drop association with any currently registered loop movers.
	void deregister_nested_loop_movers();

	/// @brief Return a reference to the list of nested loop movers.
	vector1<LoopMoverOP> const & get_nested_loop_movers() const;

private:
	bool trust_fold_tree_;
	bool manage_score_function_;
	bool was_successful_;

	Loops loops_;
	core::scoring::ScoreFunctionOP score_function_;
	vector1<LoopMoverOP> nested_movers_;
};

}
}

#endif
