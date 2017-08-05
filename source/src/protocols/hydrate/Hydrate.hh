// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hydrate/Hydrate.cc
/// @brief The Hydrate Protocol
/// @detailed
/// @author Joaquin Ambia, Jason K. Lai

#ifndef INCLUDED_protocols_hydrate_Hydrate_HH
#define INCLUDED_protocols_hydrate_Hydrate_HH

// Protocol headers
#include <protocols/hydrate/Hydrate.fwd.hh>
#include <protocols/moves/Mover.hh>

// Core and Basic headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>
#include <string>

namespace protocols {
namespace hydrate {

/// @brief
class Hydrate: public protocols::moves::Mover {
public:
	// Default constructor
	Hydrate();

	// Constructor that allows to specity a score function and the protein flexibility (region to pack & minimize )
	Hydrate(
		core::scoring::ScoreFunctionOP scorefxn,
		std::string protein_flexibility = "all"
	);

	// Copy constructor
	Hydrate(Hydrate const & hyd);

	// Assignment operator
	Hydrate & operator=(Hydrate const & hyd);

	// Destructor
	virtual ~Hydrate();

	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief  Apply the corresponding move to <input_pose>.
	virtual void apply(core::pose::Pose & input_pose);


private:
	// Initialize data members.
	void init();

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(Hydrate hyd_to, Hydrate hyd_from);

	core::scoring::ScoreFunctionOP score_fxn_; // Score function used in the protocol
	core::pack::task::TaskFactoryOP main_task_factory_; // Main task factory for the protocol
	std::string protein_flexibility_;     // Option to determine the regions of the protein to pack & minimize
	bool hydrate_all_;     // Triggers a function that searches for possible water places over the entire protein
	core::Real partial_hydrate_dew_; // Fraction of water molecules that will be packed during the first step using rotamers
	// with two optimized hb (double edge water, dew); defaults to 0.75
	bool short_residence_time_mode_; // Triggers a different algorithm to add water
	core::Real near_water_threshold_;  // Distance between a wat molecule and a residue to consider it flexible; yumeng
	bool minimize_bb_where_packing_;  // The minimizer will also minimize the backbone in the regions that were packed; yumeng
};  // class Hydrate

}  // namespace hydrate
}  // namespace protocols

#endif  // INCLUDED_protocols_hydrate_Hydrate_HH
