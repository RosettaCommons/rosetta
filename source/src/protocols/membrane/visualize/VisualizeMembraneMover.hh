// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file     protocols/membrane/VisualizeMembraneMover.hh
/// @brief      Visualize Membrane Planes by many atoms
/// @details    This does not represent the membrane planes as planes but rather
///    as a large number of additional HETATOMs in the PDB file.
///    IF YOU USE PYMOL, IT'S BETTER TO USE THE PYMOLMOVER INSTEAD!
///    If you use Chimera or alternate methods for visualization, it
///    it is still useful.
///    Last Modified: 6/19/14
/// @author  Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_hh
#define INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_hh

// Unit Headers
#include <protocols/membrane/visualize/VisualizeMembraneMover.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/conformation/Residue.fwd.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace membrane {
namespace visualize {

/// @brief Add membrane planes to the pose represented by
///   2 layers of MEM virtual residues
class VisualizeMembraneMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief   Defualt Constructor
	/// @details  Construct membrane residues with spacing = 5,
	///           width = 100, and the pose membrane and center
	VisualizeMembraneMover();

	/// @brief    Construct with User specified spacing & width
	/// @details  Construct membranes with a given spacing and
	///     width in angstroms
	VisualizeMembraneMover(
		core::Real spacing,
		core::Real width,
		core::Real thicnkess
	);

	/// @brief Copy Constructor
	/// @details Creates a deep copy of the visualize membrane mover class
	VisualizeMembraneMover( VisualizeMembraneMover const & src );

	/// @brief Assignment Operator
	/// @details Overloads "=" assignemnt for deep copying
	VisualizeMembraneMover &
	operator=( VisualizeMembraneMover const & src );

	/// @brief Destructor
	~VisualizeMembraneMover();

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

	/// @brief    Apply Visualize Transform "Move"
	/// @details  Adds a series of virtiaul residues to the pose given a
	///     spacing and width specified at construction time
	virtual void apply( core::pose::Pose & pose );

	/// @brief   Return the name of this mover
	virtual std::string get_name() const;

private:

	//////////////////////
	/// Helper Methods ///
	//////////////////////

	/// @brief Register Options with JD2
	void register_options();

	/// @brief Initialize Options from the Command Line
	/// @details Options allowed are vrt spacing and plane width
	void init_from_cmd();

	/// @brief Create a Membrane Residue
	/// @details Given a centered position and residue typeset, return
	/// a ResidueOP with the xyz coordinate pos, type MEM, from typeset given
	core::conformation::ResidueOP
	create_membrane_virtual( core::Vector pos, bool fullatom );

private:

	// Spacing and width of VRTs defining planes
	core::Real spacing_;
	core::Real width_;

	// Set membrane thicnkess to visualize
	core::Real thickness_;

};

} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_hh
