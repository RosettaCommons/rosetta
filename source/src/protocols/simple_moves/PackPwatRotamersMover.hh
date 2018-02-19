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
/// @author modified from PackRotamersMover.hh by Ryan Pavlovicz

#ifndef INCLUDED_protocols_simple_moves_PackPwatRotamersMover_hh
#define INCLUDED_protocols_simple_moves_PackPwatRotamersMover_hh

// Unit headers
#include <protocols/simple_moves/PackPwatRotamersMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

//#ifdef __clang__
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
//#endif
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/prepack_pwat_rotamers.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

#include <numeric/xyzVector.hh>

namespace protocols {
namespace simple_moves {

typedef utility::vector1< utility::vector1< core::Vector > > SetofSets;

/// @brief A protocols::moves::Mover that packs the side-chains using a rotamer library
/// It uses a ScoreFunction for packing and a PackerTask,
/// or a TaskFactory that generates a PackerTask, for instructions on
/// what rotamer sets are allowed at each residue position during packing
///
/// Common Methods:
///     PackRotamersMover.apply
/// @note please derive from PackRotamersMover instead of attempting to add protocol-specific stuff here!
class PackPwatRotamersMover : public protocols::moves::Mover {
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
	PackPwatRotamersMover();

	/// @brief constructor with typename
	PackPwatRotamersMover( std::string const & );

	/// @brief Constructs a PackRotamersMover with PackerTask  <task>
	/// evaluated using  <scorefxn>
	///
	/// ScoreFunction  scorefxn   /function to minimize while changine rotamers
	/// PackerTask     task       /object specifying what to design/pack
	/// Size (int)     nloop      /number of loops in the Pose (???)
	PackPwatRotamersMover(
		ScoreFunctionCOP scorefxn,
		PackerTaskCOP task = 0,
		core::Size nloop = 1
	);

	static std::string mover_name();

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~PackPwatRotamersMover();

	// copy constructor
	PackPwatRotamersMover( PackPwatRotamersMover const & other );

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

	///@brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	///@brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_score_function(
		TagCOP,
		basic::datacache::DataMap const &,
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

	core::Vector
	centroid_MW( utility::vector1< PointDwell > point_group );

	PointDwell
	centroid_MW_PointDwell( utility::vector1< PointDwell > point_group );

	PointDwell
	most_visited( utility::vector1< PointDwell > point_group );

	utility::vector1< utility::vector1< PointDwell > >
	cluster_rotset( utility::vector1< PointDwell > cutoff_set );

	Size
	find_closest( Pose const & pose, core::Vector Ocoord );

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
	core::Real pack_temp() const;
	bool cluster_results() const;
	bool limit_waters() const;
	bool lkb_pwat() const;
	bool exclude_exposed() const;
	bool use_average() const;

	static
	utility::tag::XMLSchemaComplexTypeGeneratorOP
	complex_type_generator_for_pack_pwat_rotamers_mover( utility::tag::XMLSchemaDefinition & xsd );

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

protected:
	///@brief update task to include all waters
	core::pack::task::PackerTaskOP
	update_task(
		Pose const & pose,
		core::pack::task::PackerTaskCOP & packer_task,
		bool force_include_current,
		bool only_water
	);

	///@brief update pose for new lkball/bb_pwat intersection rotamer sets
	void build_lkb_rotsets( Pose & pose, SetofSets const & new_pwat_rotsets );

	///@brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
	virtual void setup( Pose & pose );

	// need a more elegant rot_to_pack implementation than this
	virtual core::PackerEnergy run(
		Pose & pose,
		//  utility::vector0< int > rot_to_pack = utility::vector0< int >(),
		utility::vector0< int > rot_to_pack,
		utility::vector1< PointDwell > & all_rot
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

	core::Real pack_temp_;
	bool cluster_results_;
	bool limit_waters_;
	bool lkb_pwat_;
	bool exclude_exposed_;
	bool use_average_;
};

std::ostream &operator<< (std::ostream &os, PackPwatRotamersMover const &mover);

// note: it is better to create new files, instead of adding additional classes here

} // moves
} // protocols

#endif
