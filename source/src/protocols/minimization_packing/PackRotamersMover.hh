// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author ashworth (current form)

#ifndef INCLUDED_protocols_minimization_packing_PackRotamersMover_hh
#define INCLUDED_protocols_minimization_packing_PackRotamersMover_hh

// Unit headers
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>

#include <core/pack/annealer/AnnealerObserver.fwd.hh>

// Utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector0.hh>

#include <list>

#include <core/scoring/annealing/RotamerSets.fwd.hh> // AUTO IWYU For RotamerSetsCOP, RotamerSetsOP

namespace protocols {
namespace minimization_packing {

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
	/// @brief default constructor; reads the nloop_ value from the global options system.
	PackRotamersMover();

	/// @brief constructor that reads from the (possibly local) options collection object
	PackRotamersMover( utility::options::OptionCollection const & options );

	/// @brief constructor with typename; reads the nloop_ value from the global options system.
	PackRotamersMover( std::string const & );

	/// @brief Constructs a PackRotamersMover with TaskFactory  <taskfactory>
	/// evaluated using  <scorefxn>
	///
	/// ScoreFunction  scorefxn     /function to minimize while changine rotamers
	/// TaskFactory    taskfactory  /object specifying what to design/pack
	/// core::Size (int)     nloop  /number of rounds to run packing
	PackRotamersMover(
		ScoreFunctionCOP scorefxn,
		TaskFactoryCOP task_factory = nullptr,
		core::Size nloop = 1
	);

	/// @brief Constructs a PackRotamersMover with PackerTask  <task>
	/// evaluated using  <scorefxn>
	///
	/// Note: The version with the TaskFactory is likely preferred,
	/// unless you have reasons for using a fixed task.
	///
	/// ScoreFunction  scorefxn   /function to minimize while changine rotamers
	/// PackerTask     task       /object specifying what to design/pack
	/// core::Size (int)     nloop      /number of rounds to run packing
	PackRotamersMover(
		ScoreFunctionCOP scorefxn,
		PackerTaskCOP task,
		core::Size nloop = 1
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	~PackRotamersMover() override;

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
	void apply( Pose & pose ) override;

	std::string get_name() const override;

	void show(std::ostream & output=std::cout) const override;

	bool task_is_valid( Pose const & pose ) const; // should this be virtual?

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &
	) override;

	/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_score_function(
		TagCOP,
		basic::datacache::DataMap const &
	);

	/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &
	);

	void
	initialize_task_factory_with_operations(
		std::list< core::pack::task::operation::TaskOperationCOP > const &
	);

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP clone() const override;

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
	virtual void task_factory( TaskFactoryCOP tf );

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

	/// @brief Returns the RotamerSets object being used for packing
	/// @details - This is only valid after a call to setup or apply.
	RotamerSetsCOP rotamer_sets() const;

	//std::string get_name() const override { return mover_name(); }
	static std::string mover_name();

	static utility::tag::XMLSchemaComplexTypeGeneratorOP complex_type_generator_for_pack_rotamers_mover( utility::tag::XMLSchemaDefinition & xsd );
	// The above was added so that this could be called in SymPackRotamersMover
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static void list_options_read( utility::options::OptionKeyList & opts );

	void
	set_annealer_observer( core::pack::annealer::AnnealerObserverOP observer ){
		observer_ = observer;
	}

protected:
	/// @brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
	virtual void setup( Pose & pose );

	// need a more elegant rot_to_pack implementation than this
	virtual core::PackerEnergy run(
		Pose & pose,
		utility::vector0< int > rot_to_pack = utility::vector0< int >()
	) const;

	/// @brief Clean up cached pose and mover data after the fact.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual void cleanup( core::pose::Pose & pose );

	virtual void note_packertask_settings( Pose const & );

private:
	void initialize_from_options( utility::options::OptionCollection const & options );

private:
	// pointers to data that are passed in
	ScoreFunctionCOP scorefxn_;
	PackerTaskCOP task_;
	core::Size nloop_;
	TaskFactoryCOP task_factory_;

	// 'really private:' packer data, actually created and owned by this class
	RotamerSetsOP rotamer_sets_;
	AnnealableGraphBaseOP ig_;

	core::pack::annealer::AnnealerObserverOP observer_;
};

std::ostream &operator<< (std::ostream &os, PackRotamersMover const &mover);

// note: it is better to create new files, instead of adding additional classes here

} // moves
} // protocols

#endif
