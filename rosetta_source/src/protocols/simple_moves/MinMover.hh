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
/// @author
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_simple_moves_MinMover_hh
#define INCLUDED_protocols_simple_moves_MinMover_hh

// Unit headers
#include <protocols/simple_moves/MinMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

#include <protocols/filters/Filter.fwd.hh>


#include <core/kinematics/MoveMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers

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
class MinMover : public protocols::moves::Mover
{
public:

	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::optimization::MinimizerOptionsOP MinimizerOptionsOP;
	typedef core::optimization::MinimizerOptionsCOP MinimizerOptionsCOP;
	typedef core::Real Real;

public:

	// default constructor
	/// @brief Constructs a MinMover
	/// minmover = protocols::simple_moves::MinMover()
	MinMover();

	MinMover( std::string const & );

	virtual ~MinMover();

	// constructor with arguments
	MinMover(
		core::kinematics::MoveMapOP movemap_in,
		ScoreFunctionCOP scorefxn_in,
		std::string const & min_type_in,
		Real tolerance_in,
		bool use_nb_list_in,
		bool deriv_check_in = false,
		bool deriv_check_verbose_in = false
	);

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	/// @brief Called by protocols::moves::MoverFactory when constructing new protocols::moves::Movers. Takes care of the specific mover's parsing.
	virtual
	void parse_my_tag(
		TagPtr const,
		protocols::moves::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	void parse_opts(
		TagPtr const,
		protocols::moves::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	void parse_chi_and_bb( TagPtr const );

	/// @brief allow non-const access to the internal minimizer options object
	virtual MinimizerOptionsOP min_options();
	/// @brief allow const access to the internal minimizer options object
	virtual MinimizerOptionsCOP min_options() const;

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
	virtual core::kinematics::MoveMapCOP movemap();
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
	virtual void score_function( core::scoring::ScoreFunction const & scorefxn_in );
	virtual ScoreFunctionCOP score_function() const;

	virtual void min_type( std::string min_type_in );
	virtual void tolerance( Real tolerance_in );
	virtual void nb_list( bool nb_list_in );
	virtual void deriv_check( bool deriv_check_in );

	Real tolerance();
	std::string min_type();
	bool nb_list();
	bool deriv_check();


	//	void threshold( Real threshold_in ) { threshold_ = threshold_in; } // TODO: can be deleted?
	//	Real threshold() { return threshold_; } // TODO: can be deleted?

	/// @brief Minimizes the DOFs of  <pose_>  specified in the MoveMap
	/// using the ScoreFunction
	///
	/// example(s):
	///     minmover.apply(pose)
	/// See also:
	///     MinMover
	///     MinMover.movemap
	///     MinMover.score_function
	virtual void apply( core::pose::Pose & pose_ );
	virtual std::string get_name() const;

	inline void cartesian( bool newval ) { cartesian_ = newval; }
	inline bool cartesian( ) const { return cartesian_; }

private:
	// data
	core::kinematics::MoveMapOP movemap_;
	ScoreFunctionCOP scorefxn_;
	MinimizerOptionsOP min_options_;
	Real threshold_;
	bool cartesian_;
};

} // moves
} // rosetta
#endif
