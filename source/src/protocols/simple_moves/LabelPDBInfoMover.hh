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
/// @author Chris King

#ifndef INCLUDED_protocols_simple_moves_LabelPDBInfoMover_hh
#define INCLUDED_protocols_simple_moves_LabelPDBInfoMover_hh

// Unit headers
#include <protocols/simple_moves/LabelPDBInfoMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

// AUTO-REMOVED #include <core/pack/interaction_graph/InteractionGraphBase.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.fwd.hh>
//#ifdef __clang__
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
//#endif
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>


namespace protocols {
namespace simple_moves {

/// @brief A protocols::moves::Mover that adds PDBInfo REMARK labels based on task operations
///
/// Common Methods:
///     LabelPDBInfoMover.apply
/// @note please derive from LabelPDBInfoMover instead of attempting to add protocol-specific stuff here!
class LabelPDBInfoMover : public protocols::moves::Mover {
public:
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;

public:
	/// @brief default constructor
	LabelPDBInfoMover();

	/// @brief constructor with typename
	LabelPDBInfoMover( std::string const & );

	LabelPDBInfoMover(
		PackerTaskCOP task,
		std::string label
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~LabelPDBInfoMover();

	// copy constructor
	LabelPDBInfoMover( LabelPDBInfoMover const & other );

	// methods

	/// @brief Performs side-chain packing based on the input PackerTask
	/// using the input ScoreFunction
	///
	///	example(s):
	///     packmover.apply(pose)
	/// See Also:
	///     PackerTask
	///     ScoreFunction
	virtual void apply( Pose & pose );

	virtual std::string get_name() const;

	virtual void show(std::ostream & output=std::cout) const;

	bool task_is_valid( Pose const & pose ) const; // should this be virtual?

	///@brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	///@brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	///@brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP fresh_instance() const;

	///@brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP clone() const;

	// setters

	/// @brief Sets the TaskFactory to  <tf>
	///
	/// example(s):
	///     packmover.task_factory(task_design)
	/// See Also:
	///     LabelPDBInfoMover
	///     LabelPDBInfoMover.task
	void task_factory( TaskFactoryCOP tf );

	/// @brief Sets the PackerTask to  <t>
	///
	/// example(s):
	///     packmover.task(task_pack)
	/// See Also:
	///     LabelPDBInfoMover
	///     LabelPDBInfoMover.task_factory
	void task( PackerTaskCOP t );
	void label( std::string l );


	// accessors

	/// @brief Returns the PackerTask
	///
	/// example(s):
	///     packmover.task()
	/// See Also:
	///     LabelPDBInfoMover
	///     LabelPDBInfoMover.task_factory
	PackerTaskCOP task() const;

	/// @brief Returns the TaskFactory
	///
	/// example(s):
	///     packmover.task_factory()
	/// See Also:
	///     LabelPDBInfoMover
	///     LabelPDBInfoMover.task
	TaskFactoryCOP task_factory() const;
	std::string label() const;

protected:

private:
	// pointers to data that are passed in
	PackerTaskCOP task_;
	TaskFactoryCOP task_factory_;
	std::string label_;
};

std::ostream &operator<< (std::ostream &os, LabelPDBInfoMover const &mover);

} // moves
} // protocols

#endif
