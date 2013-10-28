// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//          under license.
// (c) The Rosetta software is developed by the contributing members of the
//          Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
//          Questions about this can be
// (c) addressed to University of Washington UW TechTransfer,
//          email: license@u.washington.edu.

/// @file FullatomRelaxMover.hh
/// @brief
/// @author Robin A Thottungal  (rathottungal@gmail.com)

#ifndef INCLUDED_protocols_surface_docking_FullatomRelaxMover_hh
#define INCLUDED_protocols_surface_docking_FullatomRelaxMover_hh

// Unit Headers
#include <protocols/surface_docking/FullatomRelaxMover.fwd.hh>
// Package headers
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>
// Project headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/docking/DockingHighRes.hh>
#include <core/types.hh>
#include <utility/vector1_bool.hh>


// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>
#include <map>
#include <list>

//Auto Headers
#include <sstream>



namespace protocols {
namespace surface_docking {

class FullatomRelaxMover : public moves::Mover {

public:

	FullatomRelaxMover();

	//destructor
	~FullatomRelaxMover();

	// virtual functions that get overloaded or
        // called from the inheriting classes
	void apply( core::pose::Pose & );

	virtual std::string get_name() const;

	//virtual void setup_list( core::pose::Pose & ) = 0;

	//virtual void set_angles( core::Real ) = 0;

	//virtual bool make_move( core::pose::Pose & ) = 0;

	// Undefined, commenting out to fix PyRosetta build  void angle_max(const core::Real);

	// Undefined, commenting out to fix PyRosetta build  void angle_max(const char,const core::Real);

	void set_nmoves(const core::Size);

	void setup_defaults();

	// Undefined, commenting out to fix PyRosetta build  void setup_shearTrialMover();

	// Undefined, commenting out to fix PyRosetta build  void setup_smallTrialMover();

	void setupMovers();

	void FinalizeMovers(core::pose::Pose &);

	void init_from_options();

	void set_smallmovesize(Size scale);

	void set_ljrepulsion_weight(core::Real weight_scale);

	void set_ecounter(core::Size ecount);

	void CalcSecondayStruct(core::pose::Pose & pose);


	std::string getSecondayStruct();

public:
	typedef core::Real Real;
private:

	// Options
	bool benchmark_;

	//members for smallTrialMove
	core::Real temperature_;
	core::Real kT_;
	core::Size nmoves_;
	core::Size encounter_;
	core::Size encounter_cycle_;

	std::map< char, core::Real > angle_max_;
        // for scoring
	core::scoring::ScoreFunctionOP score_high_res_;
	// for movers

	std::string smallmin_type_;
	std::string shearmin_type_;

	core::kinematics::MoveMapOP moveMapOP_;
	core::kinematics::MoveMapOP smallmoveMapOP_;
	core::kinematics::MoveMapOP shearmoveMapOP_;

	moves::MonteCarloOP  monteCarlo_;
	moves::MonteCarloOP  smallmonteCarlo_;
	moves::MonteCarloOP  shearmonteCarlo_;

	simple_moves::SmallMoverOP smallmover_;
	simple_moves::MinMoverOP smallminmover_;
	moves::SequenceMoverOP smallsequenceMover_;
	moves::TrialMoverOP small_trial_min_mover_;

	simple_moves::ShearMoverOP shearmover_;
	simple_moves::MinMoverOP shearminmover_;
	moves::SequenceMoverOP shearsequenceMover_;
	moves::TrialMoverOP shear_trial_min_mover_;


	//protocols::docking::DockingHighResOP dockingHigh_res_;
	std::string SecStruct_;
    };


} // surfaceDockingProtocols
} // protocols

#endif
