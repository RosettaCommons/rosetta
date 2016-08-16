// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Optimizes the membrane position given the high-res score function
/// @details Optimizes the membrane position given the smooth high-res score
///    function; scans the center along the normal around the initial center
///    in 0.1A steps; scans the normal in 0.2 degree steps along arches
///    over the x-axis, y-axis, xy-direction, -xy-direction; outcome is
///    deterministic
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_OptimizeMembranePositionMover_hh
#define INCLUDED_protocols_membrane_OptimizeMembranePositionMover_hh

// Unit Headers
#include <protocols/membrane/OptimizeMembranePositionMover.fwd.hh>
#include <protocols/membrane/OptimizeMembranePositionMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

class OptimizeMembranePositionMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Defaults: scorefxn = smooth2012
	OptimizeMembranePositionMover();

	/// @brief Copy Constructor
	OptimizeMembranePositionMover( OptimizeMembranePositionMover const & src );

	/// @brief Assignment Operator
	OptimizeMembranePositionMover & operator = ( OptimizeMembranePositionMover const & src );

	/// @brief Destructor
	virtual ~OptimizeMembranePositionMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (OptimizeMembranePositionMover)
	virtual std::string get_name() const;

	/// @brief Flip the downstream partner in the membrane
	virtual void apply( core::pose::Pose & pose );

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	/// @brief Set default values
	void set_defaults();

	/// @brief optimize membrane center
	void optimize_membrane_center( core::pose::Pose & pose );

	/// @brief Optimize membrane normal
	void optimize_membrane_normal( core::pose::Pose & pose );

private: // data

	/// @brief Original foldtree before optimization
	/// @details Will be reset to original after the optimization
	core::kinematics::FoldTree ft_;

	/// @brief Scorefunction
	core::scoring::ScoreFunctionOP sfxn_;
	core::Real score_best_;

	/// @brief center search
	core::Real starting_z_;
	core::Real best_z_;
	core::Real stepsize_z_;
	core::Vector best_center_;

	/// @brief Normal search
	core::Real stepsize_angle_;
	core::Vector best_normal_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_OptimizeMembranePositionMover_hh
