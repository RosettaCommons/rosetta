// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/visualize/ShowMembranePlanesMover.hh
///
/// @brief 		Add Anchor Residues for Membrane Planes
/// @details    Add 6 virtual membrane residues to the pose as an additional
///				chain to the protein. Anchor residues are attached by jump to the membrane center
///				and are pointing downstream to allow movement of the planes. This plug in
///				works directly with the PyMol Mover
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 6/27/14

#ifndef INCLUDED_protocols_membrane_visualize_ShowMembranePlanesMover_hh
#define INCLUDED_protocols_membrane_visualize_ShowMembranePlanesMover_hh

// Unit Headers
#include <protocols/membrane/visualize/ShowMembranePlanesMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/conformation/Residue.fwd.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace membrane {
namespace visualize {

using namespace core;
using namespace core::conformation;
using namespace protocols::moves;

/// @brief Add Anchor Residues for Representing Membrane Planes
/// Prerequisite mover for using the PyMol Observer/Viewer in Rosetta
class ShowMembranePlanesMover : public Mover {
	
public:
	
	/////////////////////
	//// Constructors ///
	/////////////////////
	
	/// @brief Default Constructor
	/// @details Construct a Show membrane planes mover - uses the membrane framework
	/// normal and center parameters to compute 6 residues defining an upper/lower plane
	/// as the membrane planes
	ShowMembranePlanesMover();
	
	/// @brief Constructor
	/// @details Construct a Show membrane planes mover - uses the membrane framework
	/// normal and center parameters to compute 6 residues defining an upper/lower plane
	/// as the membrane planes. Can specify membrane thicnkess and/or raidus of gyration to define
	/// the radii of the planes (triangle shaped)
	ShowMembranePlanesMover(
							core::Real thicnkess,
							core::Size npoints
							);
	
	/// @brief Copy Constructor
	/// @details Create a deep copy of the show membrane planes mover
	ShowMembranePlanesMover( ShowMembranePlanesMover const & src );
	
	/// @brief Assignemnt Operator
	/// @details Create a deep copy of the Show membrane planes mover while
	/// overloading the assignment operator
	ShowMembranePlanesMover &
	operator=( ShowMembranePlanesMover const & src );
	
	/// @brief Destructor - Get rid of this class
	~ShowMembranePlanesMover();
	
	//////////////////////
	//// Mover Methods ///
	//////////////////////
	
	/// @brief Get the name of this mover (ShowMembranePlanesMover)
	virtual std::string get_name() const;
	
	/// @brief Add Anchor Residues to Pose
	/// @details Add 6 anchor residues to pose defining the bounds of the membrane planes.
	/// Width of membrane planes are determined by 2x RG of the pose.
	virtual void apply( Pose & pose );
	
	////////////////////////////////
	//// Rosetta Scripts Methods ///
	////////////////////////////////
	
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
	ResidueOP
	create_membrane_virtual( Vector pos, bool fullatom );
	
	/// @brief Select Points to Define the Upper and Lower Planes
	utility::vector1< Vector >
	select_plane_points(
		core::Vector center,
		core::Vector normal,
		core::Size n_points,
		core::Real radius
		);

private: 
	
	// Membrane thickness
	core::Real thickness_;
	
	// Number of points to define the plane
	core::Size npoints_;
	
	
};

} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_ShowMembranePlanesMover_hh
