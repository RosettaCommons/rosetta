// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/MinMover.hh
/// @brief  class definition for MinMover
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_simple_moves_MinMover_HH
#define INCLUDED_protocols_simple_moves_MinMover_HH

// Unit headers
#include <protocols/simple_moves/MinMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoveMapMover.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace simple_moves {

///////////////////////////////////////////////////////////////////////////////
/// @brief A protocols::moves::Mover that minimizes a Pose to a local energy minimum by
/// performing energy minimization of a ScoreFunction over the allowable
/// degrees of freedom, defined by a MoveMap. The minimization type,
/// minimization tolerance, and various other options can be also be set.
///
/// Common Methods:
///     MinMover.apply
///     MinMover.movemap
///     MinMover.score_function
class MinMover : public protocols::moves::MoveMapMover {
	typedef protocols::moves::MoveMapMover Parent;
public:

	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::optimization::MinimizerOptionsOP MinimizerOptionsOP;
	typedef core::optimization::MinimizerOptionsCOP MinimizerOptionsCOP;
	typedef core::Real Real;

private:
	// This type needs both DOF_Type and Torsion type to handle both
	// backbone torsion angles and bond angles and bond lengths.
	typedef
		std::map<
		std::pair< core::id::DOF_Type, core::id::TorsionType >,
		core::pack::task::TaskFactoryOP >
		DOF_TaskMap;

public:
	// default constructor
	/// @brief Constructs a MinMover
	/// minmover = protocols::simple_moves::MinMover()
	MinMover();

	MinMover( std::string const & );

	~MinMover() override;

	/// @brief constructor with arguments
	///
	/// @details
	/// min_type_in indicates the type of minimizer to use.  This is different for cartesion.
	///  Current recommendation (Jan 2016) is "dfpmin_armijo_nonmonotone"
	///
	/// tolerance_in indicates the tolerance on the minimizer.  Smaller tolerance will result
	///  in lower energies, but increased runtime.  Generally one would use .01
	///
	/// use_nb_list_in The neighbor list is based on the move map.  It includes any
	///  atoms that can be moved by the minimizer plus their neighbors.  This list
	///  is not updated during minimization.  All scores for atoms and atom pairs
	///  outside the neighbor list are held fixed.  All scores for atoms and atom
	///  pairs within the list are not cached (usually they would be), since it
	///  assumed that they will be changing rapidly.  These optimizations are
	///  effective when a large number of small moves are being made.  You may
	///  prefer to enable this option when minimizing in fullatom mode, and to
	///  disable it in centroid mode.
	///
	/// @see core::scoring::AtomNeighbor
	///
	MinMover(
		core::kinematics::MoveMapOP movemap_in,
		ScoreFunctionCOP scorefxn_in,
		std::string const & min_type_in,
		Real tolerance_in,
		bool use_nb_list_in,
		bool deriv_check_in = false,
		bool deriv_check_verbose_in = false
	);

	/// @brief Minimizes the DOFs of  <pose_>  specified in the MoveMap
	/// using the ScoreFunction
	///
	/// example(s):
	///     minmover.apply(pose)
	/// See also:
	///     MinMover
	///     MinMover.movemap
	///     MinMover.score_function
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return mover_name(); }
	static std::string mover_name();
	void show(std::ostream & output=std::cout) const override;

	inline void cartesian( bool newval ) { cartesian_ = newval; }
	inline bool cartesian( ) const { return cartesian_; }

	inline void abs_score_convergence_threshold( Real newval ) { abs_score_convergence_threshold_ = newval; }
	inline Real abs_score_convergence_threshold( ) const { return abs_score_convergence_threshold_; }

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Called by protocols::moves::MoverFactory when constructing new protocols::moves::Movers. Takes care of the specific mover's parsing.

	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) override;

	void parse_opts(
		TagCOP,
		basic::datacache::DataMap & data_map,
		Filters_map const &,
		protocols::moves::Movers_map const & );

	void parse_movemap_factory( TagCOP, basic::datacache::DataMap & );

	void parse_dof_tasks(
		TagCOP tag,
		basic::datacache::DataMap & data);

	void parse_dof_task_type(
		std::string const & tag_name,
		core::id::DOF_Type dof_type,
		core::id::TorsionType torsion_type,
		TagCOP tag,
		basic::datacache::DataMap & data);

	/// @brief allow non-const access to the internal minimizer options object
	virtual MinimizerOptionsOP min_options();

	/// @brief allow const access to the internal minimizer options object
	virtual MinimizerOptionsCOP min_options() const;

	/// @brief directly set the internal minimizer options object
	virtual void min_options(MinimizerOptionsOP min_options);

	/// @brief Sets the MoveMap to  <movemap_in>
	/// determines which DOF to minimize
	///
	/// example(s):
	///     minmover.movemap(movemap1)
	/// See also:
	///     MinMover
	///     MinMover.apply
	///     MinMover.score_function
	///     MoveMap
	virtual void movemap( core::kinematics::MoveMapCOP movemap_in );
	void set_movemap( core::kinematics::MoveMapCOP movemap_in ) override;

	/// @brief Sets the MoveMapFactory
	/// The MoveMapFactory will be used to construct a MoveMap if an explicit one isn't set
	virtual void movemap_factory( core::select::movemap::MoveMapFactoryCOP mmf );

	/// @brief Get the movemap for the pose
	/// @details This is the preferred version to call.
	core::kinematics::MoveMapCOP movemap( core::pose::Pose const & pose ) const override;

	/// @brief Get the explicitly set movemap object, if present, nullptr if not
	/// Will not create a default or look at the MoveMapFactory
	core::kinematics::MoveMapCOP explicitly_set_movemap() const;

	/// @brief Sets the ScoreFunction to  <scorefxn_in>
	/// determines which function to minimize
	///
	/// example(s):
	///     minmover.score_function(scorefxn)
	/// See also:
	///     MinMover
	///     MinMover.apply
	///     MinMover.movemap
	///     ScoreFunction
	virtual void score_function( ScoreFunctionCOP scorefxn_in );
	//virtual void score_function( core::scoring::ScoreFunction const & scorefxn_in );
	virtual ScoreFunctionCOP score_function() const;

	virtual void min_type( std::string min_type_in );
	virtual void tolerance( Real tolerance_in );
	/// @copydoc core::optimization::MinimizerOptions::use_nblist(bool)
	virtual void nb_list( bool nb_list_in );
	virtual void deriv_check( bool deriv_check_in );

	virtual void max_iter( Size max_iter_in );

	Real tolerance() const;
	std::string min_type() const;
	/// @copydoc core::optimization::MinimizerOptions::use_nblist()
	bool nb_list() const;
	bool deriv_check() const;
	bool omega() const{ return omega_; }
	void omega( bool const b ){ omega_ = b; }
	Real score_diff_after_minimization() const;
	Real abs_score_diff_after_minimization() const;

	static utility::tag::XMLSchemaComplexTypeGeneratorOP complex_type_generator_for_min_mover( utility::tag::XMLSchemaDefinition & xsd );
	// The above was added so that this could be called in SymMinMover
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

protected:

	/// @brief for use with RosettaScripts current method of using
	///TaskOperations to specify residue sets.
	void
	apply_dof_tasks_to_movemap(
		core::pose::Pose const & pose,
		core::kinematics::MoveMap & movemap) const;

	/// @brief run the minimizer (cartesian or torsion) based on set options, with the final movemap setup. Will call inner_run_minimizer multiple times if
	/// abs_score_convergence_threshold is set
	void
	minimize( core::pose::Pose & pose, core::kinematics::MoveMap & active_movemap );

private:
	/// @brief actual call to the minimizer, in helper function so it can be
	/// run multiple times by minimize()
	void
	inner_run_minimizer( core::pose::Pose & pose, core::kinematics::MoveMap & active_movemap );

	// data

	// The MoveMapFactory will be overridden by an explicitly set movemap_.
	core::select::movemap::MoveMapFactoryCOP movemap_factory_;
	core::kinematics::MoveMapOP movemap_;
	bool omega_; //dflt true ; minimize omega?
	ScoreFunctionCOP scorefxn_;
	MinimizerOptionsOP min_options_;
	Real threshold_;
	Real abs_score_convergence_threshold_;
	Real score_before_minimization_;
	Real score_after_minimization_;
	bool cartesian_;

	/// @details Until ResidueSubsetOperations are implemented,
	///RosettaScripts uses TaskOperations as a generic way to specify
	///sets of residues, so it is necessary to have a method of
	///controling a movemap from a set of task operations. This is not
	///suppose to be a general interface for MinMover, so only expose it
	///through the parse_my_tag. Another small complication is that
	///TaskOperations must be defined in the context of a pose, so hang
	///on to it here until the apply. In case MinMover is applied
	///multiple times, to avoid accumulating state, reset the move map
	///to before the task operations were applied at the end of apply.
	DOF_TaskMap dof_tasks_;
};

std::ostream &operator<< (std::ostream &os, MinMover const &mover);

}  // simple_moves
}  // protocols

#endif  // INCLUDED_protocols_simple_moves_MinMover_HH
