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
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author ashworth (current form)

#ifndef INCLUDED_protocols_simple_moves_PackRotamersMover_hh
#define INCLUDED_protocols_simple_moves_PackRotamersMover_hh

// Unit headers
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

#include <core/pack/task/PackerTask.fwd.hh>
//#ifdef __clang__
//#endif
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>


namespace protocols {
namespace simple_moves {

/// @brief A protocols::moves::Mover that packs the side-chains using a rotamer library
/// It uses a ScoreFunction for packing and a PackerTask,
/// or a TaskFactory that generates a PackerTask, for instructions on
/// what rotamer sets are allowed at each residue position during packing
///
/// Common Methods:
///     PackRotamersMover.apply
/// @note please derive from PackRotamersMover instead of attempting to add protocol-specific stuff here!
class PackRotamersMover : public protocols::moves::Mover {
public:
	typedef core::pack::interaction_graph::AnnealableGraphBaseOP AnnealableGraphBaseOP;
	typedef core::pack::interaction_graph::AnnealableGraphBaseCOP AnnealableGraphBaseCOP;
	typedef core::pack::rotamer_set::RotamerSetsOP RotamerSetsOP;
	typedef core::pack::rotamer_set::RotamerSetsCOP RotamerSetsCOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;

public:
	/// @brief default constructor
	PackRotamersMover();

	/// @brief constructor with typename
	PackRotamersMover( std::string const & );

	/// @brief Constructs a PackRotamersMover with PackerTask  <task>
	/// evaluated using  <scorefxn>
	///
	/// ScoreFunction  scorefxn   /function to minimize while changine rotamers
	/// PackerTask     task       /object specifying what to design/pack
	/// Size (int)     nloop      /number of rounds to run packing
	PackRotamersMover(
		ScoreFunctionCOP scorefxn,
		PackerTaskCOP task = 0,
		core::Size nloop = 1
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~PackRotamersMover();

	// copy constructor
	PackRotamersMover( PackRotamersMover const & other );

	// methods

	/// @brief Performs side-chain packing based on the input PackerTask
	/// using the input ScoreFunction
	///
	/// example(s):
	///     packmover.apply(pose)
	/// See Also:
	///     PackerTask
	///     ScoreFunction
	virtual void apply( Pose & pose );

	virtual std::string get_name() const;

	virtual void show(std::ostream & output=std::cout) const;

	bool task_is_valid( Pose const & pose ) const; // should this be virtual?

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_score_function(
		TagCOP,
		basic::datacache::DataMap const &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP clone() const;

	// setters

	/// @brief Sets the ScoreFunction to  <sf>
	///
	/// example(s):
	///     packmover.score_function(scorefxn)
	/// See Also:
	///     PackRotamersMover
	///     PackRotamersMover.task
	void score_function( ScoreFunctionCOP sf );

	/// @brief Sets the TaskFactory to  <tf>
	///
	/// example(s):
	///     packmover.task_factory(task_design)
	/// See Also:
	///     PackRotamersMover
	///     PackRotamersMover.task
	void task_factory( TaskFactoryCOP tf );

	/// @brief Sets the PackerTask to  <t>
	///
	/// example(s):
	///     packmover.task(task_pack)
	/// See Also:
	///     PackRotamersMover
	///     PackRotamersMover.task_factory
	void task( PackerTaskCOP t );
	void nloop( core::Size nloop_in );


	// accessors

	/// @brief Returns the ScoreFunction
	///
	/// example(s):
	///     packmover.score_function()
	/// See Also:
	///     PackRotamersMover
	///     PackRotamersMover.task
	ScoreFunctionCOP score_function() const;

	/// @brief Returns the PackerTask
	///
	/// example(s):
	///     packmover.task()
	/// See Also:
	///     PackRotamersMover
	///     PackRotamersMover.task_factory
	PackerTaskCOP task() const;

	core::Size nloop() const { return nloop_; }

	/// @brief Returns the TaskFactory
	///
	/// example(s):
	///     packmover.task_factory()
	/// See Also:
	///     PackRotamersMover
	///     PackRotamersMover.task
	TaskFactoryCOP task_factory() const;
	RotamerSetsCOP rotamer_sets() const;
	AnnealableGraphBaseCOP ig() const;

protected:
	/// @brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
	virtual void setup( Pose & pose );

	// need a more elegant rot_to_pack implementation than this
	virtual core::PackerEnergy run(
		Pose & pose,
		utility::vector0< int > rot_to_pack = utility::vector0< int >()
	) const;
	virtual void note_packertask_settings( Pose const & );

private:
	// pointers to data that are passed in
	ScoreFunctionCOP scorefxn_;
	PackerTaskCOP task_;
	core::Size nloop_;
	TaskFactoryCOP task_factory_;

	// 'really private:' packer data, actually created and owned by this class
	RotamerSetsOP rotamer_sets_;
	AnnealableGraphBaseOP ig_;
};

std::ostream &operator<< (std::ostream &os, PackRotamersMover const &mover);

// note: it is better to create new files, instead of adding additional classes here

} // moves
} // protocols

#endif
