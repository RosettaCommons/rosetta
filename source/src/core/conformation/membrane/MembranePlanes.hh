// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/conformation/MembranePlanes.hh
///
/// @brief 		Specification for Membrane Planes Residue Definitions
///	@details	When using the membrane code, users can optionally view the membrane planes
///				defined by the center/normal positions. This object will store the positions of anchoring residues
///				which will be used to draw CGO planes in PyMOl, representing the membrane planes.
///				This data should not be used in dynamic simulations, it's only purpose is for visualization.
///				Last Modified: 7/23/14
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembranePlanes_hh
#define INCLUDED_core_conformation_membrane_MembranePlanes_hh

// Unit Headers
#include <core/conformation/membrane/MembranePlanes.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <iostream>

namespace core {
namespace conformation {
namespace membrane {

/// @brief Store location of top and bottom marker residues for visualizing
/// location of the membrane planes in a simulation. Accessible via membrane info
class MembranePlanes : public utility::pointer::ReferenceCount {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////
	
	/// @brief Custom Constructor
	/// @details Create a MembranePlanes object detailing the positions of
	/// the upper and lower membrane planes
	MembranePlanes(
				   utility::vector1< Size > top_points,
				   utility::vector1< Size > bottom_points
				   );
	
	/// @brief Copy Constructor
	/// @details Make a deep copy of this membrane planes object
	MembranePlanes( MembranePlanes const & src );
	
	/// @brief Assignment Operator
	/// @details Make a deep copy of this obejct overriding "="
	MembranePlanes &
	operator=( MembranePlanes const & src );
	
	/// @brief Destructor
	~MembranePlanes();
	
	/// @brief Show Membrane Planes Info Data
	virtual void show( std::ostream & output=std::cout ) const;
	
	///////////////////////////
	/// Data Access Methods ///
	///////////////////////////
	
	/// @brief Access point residues defining top membrane plane
	utility::vector1< Size >
	top_points();
	
	/// @brief Access point residues defining bottom membrane planes
	utility::vector1< Size >
	bottom_points();
	
private: // methods
	
	/// @brief Default Constructor
	MembranePlanes();
	
private:
	
	// Data
	utility::vector1< Size > top_points_;
	utility::vector1< Size > bottom_points_;
	
};

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembranePlanes_hh
