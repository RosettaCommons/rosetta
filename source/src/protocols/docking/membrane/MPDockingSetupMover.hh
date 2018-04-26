// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/MPDockingSetupMover.fwd.hh
/// @brief      Reads in 2 poses and 2 spanfiles, concatenates them, and
///    prints them out
///    CURRENTLY ONLY WORKS FOR 2 POSES!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (10/16/14)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingSetupMover_hh
#define INCLUDED_protocols_docking_membrane_MPDockingSetupMover_hh

// Unit Headers
#include <protocols/docking/membrane/MPDockingSetupMover.fwd.hh>
#include <protocols/moves/Mover.hh>

//TAKE OUT!!!
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.fwd.hh>

// Project Headers
//#include <protocols/membrane/AddMembraneMover.fwd.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

namespace protocols {
namespace docking {
namespace membrane {

/// @brief Setup for MPDocking
/// @details Reads in 2 poses and 2 spanfiles, concatenates the poses and
///    spanfiles and prints them out
///    CURRENTLY ONLY WORKS FOR 2 POSES!!!
class MPDockingSetupMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	MPDockingSetupMover();

	/// @brief Copy Constructor
	MPDockingSetupMover( MPDockingSetupMover const & src );

	/// @brief Destructor
	virtual ~MPDockingSetupMover();

public: // methods

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Parse Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag
	);

	/////////////////////
	/// Get Methods   ///
	/////////////////////

	/// @brief get optimize1
	bool get_optimize1() const;

	/// @brief get optimize1
	bool get_optimize2() const;

	/// @brief get pose vector
	utility::vector1< core::pose::PoseOP > get_poses() const;

	/// @brief get spanfile vector
	utility::vector1< std::string > get_spanfiles() const;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (MPDockingSetupMover)

	/// @brief Reads in 2 poses and 2 spanfiles, concatenates the poses and
	///   spanfiles and prints them out
	void apply( Pose & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // methods

	// read poses
	void read_poses();

	// read spanfiles
	void read_spanfiles();

	// initialize from commandline
	void init_from_cmd();

	// transform pose into membrane
	void transform_pose_into_membrane( core::pose::Pose & pose, core::Vector center, core::Vector normal, std::string spanfile, core::Size partner );

private: // data

	// vector of poses and spanfiles that are used as input to create final pose and spanfile
	utility::vector1< core::pose::PoseOP > poses_;
	utility::vector1< std::string > spanfiles_;

	// should the position of partner 1 or 2 be optimized in the membrane?
	bool optimize1_;
	bool optimize2_;

	// Tell the compiler that we are not hiding the base
	// function with the parse_my_tag written above
	using protocols::moves::Mover::parse_my_tag;

};

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingSetupMover_hh
