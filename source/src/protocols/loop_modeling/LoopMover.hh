// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopMover_HH
#define INCLUDED_protocols_loop_modeling_LoopMover_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.fwd.hh>
#include <protocols/loop_modeling/LoopModelerTests.fwd.hh> //for friendship

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/datacache/HierarchicalDataMap.hh>
#include <boost/utility.hpp>

// RosettaScripts headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

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

/// @brief Key names for data shared between loop movers.
struct ToolboxKeys {
	static const std::string LOOPS;
	static const std::string SCOREFXN;
	static const std::string TASK_FACTORY;
};

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

class LoopMover : public moves::Mover {

public:

	/// @brief Default constructor.
	LoopMover();

	/// @brief Configure from a RosettaScripts tag.
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose);

	/// @brief Sample the pose in the regions specified by get_loops().
	/// @details The parent class apply() method automatically sets up a fold
	/// tree (if necessary) and keeps track of whether or not the move succeeded.
	/// Child classes should reimplement do_apply() instead of this method.
	void apply(Pose & pose);

	/// @brief Return the name of this mover.
	virtual string get_name() const { return "LoopMover"; }

public:

	/// @brief Add the names of all the algorithms invoked by this loop mover to
	/// the given list.  Indentation is used to represent hierarchy.
	virtual void get_children_names(
		utility::vector1<std::string> & names, std::string indent="") const;

	/// @brief Return true if the previous move was successful.
	bool was_successful() const;

	/// @brief Return the loops to be sampled on the next call to apply().
	LoopsOP get_loops();

	/// @brief Return the loops to be sampled on the next call to apply().
	LoopsCOP get_loops() const;

	/// @brief Return the specified loop.
	Loop const & get_loop(Size index) const;

	/// @brief Set the loops to be sampled on the next call to apply().
	void set_loops(LoopsOP loops);

	/// @brief Set the loops to be sampled on the next call to apply().
	void set_loops(Loops const & loops);

	/// @brief Set the loop to be sampled on the next call to apply().
	void set_loop(Loop const & loop);

	/// @brief Request a tool from this mover or any of its parents.
	/// @details If the tool cannot be found, an exception will be thrown.
	template <typename ToolTypeOP>
	ToolTypeOP get_tool(std::string key) const {
		return toolbox_->get<ToolTypeOP>(key, "");
	}

	/// @brief Request a tool from this mover or any of its parents.
	/// @details If the tool isn't be found, the given fallback will be returned.
	template <typename ToolTypeOP>
	ToolTypeOP get_tool(std::string key, ToolTypeOP fallback) const {
		return toolbox_->get<ToolTypeOP>(key, "", fallback);
	}

	/// @brief Provide a tool for this mover or any of its children to use.
	template <typename ToolTypeOP>
	ToolTypeOP set_tool(std::string key, ToolTypeOP value) {
		return toolbox_->set<ToolTypeOP>(key, "", value);
	}

public:

	/// @brief Return an enum representing the kind of fold tree that is
	/// compatible with this mover.
	/// @details The FoldTreeRequest enum values can be combined using the
	/// bitwise logical operators.  For example, you can request either the
	/// standard fold tree or a simple fold tree with `FTR_LOOPS_WITH_CUTS |
	/// FTR_SIMPLE_TREE.`
	virtual FoldTreeRequest request_fold_tree() const;

	/// @brief Promise that the calling code will setup a fold tree compatible
	/// with request_fold_tree().  If this method is not called, this mover will
	/// setup a fold tree on its own every time apply() is called.
	void trust_fold_tree();

	/// @brief Setup the given pose with a fold tree that is compatible with the
	/// given loops and requests.
	static void setup_fold_tree(
		Pose & pose, LoopsCOP loops, FoldTreeRequest request);

protected:

	/// @brief Perform the loop sampling move.  This method can be overwritten in
	/// child classes.
	virtual bool do_apply(Pose & pose);

	/// @brief Perform the loop sampling move.  This method can be overwritten in
	/// child classes.
	virtual bool do_apply(Pose & pose, Loop const & loop);

	/// @brief Add a child to this mover.
	/// @details The purpose of calling add_child() is to make it easy to use the
	/// child from within the do_apply() method of the parent.  Child movers are
	/// always configured with the same loops and the same fold tree as their
	/// parents.  Their toolboxes are linked to their parent's as well, so that
	/// the children's tools will default to the parent's tools.  Fold tree
	/// requests made by child movers will be taken into account by the parent.
	template <typename ChildSubclassOP>
	ChildSubclassOP add_child(ChildSubclassOP child) {
		using utility::pointer::dynamic_pointer_cast;

		// Refuse to add a child that already has a parent.
		if ( child->parent_name_ != "" ) {
			std::stringstream message;
			message << "Cannot add child to " << get_name() << " because child (" <<
				child->get_name() << ") already has a parent (" <<
				child->parent_name_ << ")." << std::endl << std::endl <<
				"If you want the same child to have two parents, you have to copy the "
				"child and give a copy to each parent.  If you want to switch parents, "
				"you have to call " << child->parent_name_ << ".remove() first." << std::endl;
			utility_exit_with_message(message.str());
		}

		child->parent_name_ = get_name();
		child->toolbox_->set_parent(toolbox_);
		children_.push_back(child);
		child->trust_fold_tree();
		return child;
	}

	/// @brief Remove a child from this mover.
	/// @see LoopMover::add_child()
	/// @see LoopMover::clear_children()
	void remove_child(LoopMoverOP child);

	/// @brief Remove all children from this mover.
	/// @see LoopMover::add_child()
	/// @see LoopMover::remove_child()
	void clear_children();

	/// @brief Return all of this mover's children.
	/// @see LoopMover::add_child()
	/// @see LoopMover::count_children()
	LoopMoverOPs get_children() const;

	/// @brief Return the number of children this mover has.
	/// @see LoopMover::add_child()
	/// @see LoopMover::get_children()
	Size count_children() const;

	friend class ::LoopModelerTests;

private:
	LoopMoverOPs children_;
	basic::datacache::HierarchicalDataMapOP toolbox_;
	string parent_name_;
	bool trust_fold_tree_;
	bool was_successful_;

};

}
}

#endif
