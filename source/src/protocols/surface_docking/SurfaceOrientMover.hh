// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SurfaceOrientMover.hh
/// @brief <add a description of the class>
/// @author Robin A Thottungal (raugust1@jhu.edu)
/// @author Michael Pacella (mpacella88@gmail.com)

#ifndef INCLUDED_protocols_surface_docking_SurfaceOrientMover_hh
#define INCLUDED_protocols_surface_docking_SurfaceOrientMover_hh

// Unit Headers
#include <protocols/surface_docking/SurfaceOrientMover.fwd.hh>

// Package headers
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// C++ Headers
#include <string>


namespace protocols {
namespace surface_docking {

class SurfaceOrientMover : public moves::Mover {

public:
	//Standard methods
	/// @brief Default constructor
	SurfaceOrientMover();

	//destructor
	~SurfaceOrientMover();

	void apply( core::pose::Pose & );

	virtual std::string get_name() const;

	void set_surface_parameters( protocols::surface_docking::SurfaceParametersOP surface_parameters);


private:
	//methods
	core::Vector calculate_recenter_vector( core::Vector const & total_displacement );

	//data
	protocols::surface_docking::SurfaceParametersOP surface_parameters_;

};


} // surfaceDockingProtocols
} // protocols

#endif
