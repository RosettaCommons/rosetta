// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/symmetry/SetupForMembraneSymmetry.hh
///
/// @brief		Setup a Symmetric Membrane Protein Using the Membrane Framework
/// @details	The setup for membrane symmetry class first adds a membrane residue
///				to the asymmetric unit of a protein, creates a symmetric complex, and
///				then adds the remainder of the membrane framework to capture the entire 
///				symmetric system. This should work for both relax and docking protocols
///				currently in Rosetta. 
///
///				Last Modified: 10/3/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetry_hh
#define INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetry_hh

// Unit Headers
#include <protocols/membrane/symmetry/SetupForMembraneSymmetry.fwd.hh> 
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh> 
#include <core/types.hh> 

// Utiility Headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace membrane {
namespace symmetry {

using namespace protocols::moves;
using namespace core::pose;

/// @brief Setup For Membrane Symmetry Mover
class SetupForMembraneSymmetry : public protocols::moves::Mover {

public: // constructors

	/// @brief Default Constructor for Setup for MP Symm Mover
	SetupForMembraneSymmetry();
	
	/// @brief Custom Constructor for Setup for MP Symm Mover
	SetupForMembraneSymmetry( std::string symmdef_file );

	/// @brief Create a Deep Copy of this setup for membrane symmetry mover
	SetupForMembraneSymmetry( SetupForMembraneSymmetry const & src );

	/// @brief Default Destructor
	~SetupForMembraneSymmetry();

public: // rosetta scripts methods

	/// @brief 
	MoverOP
	clone() const; 

	/// @brief 
	MoverOP
	fresh_instance() const; 

	/// @brief 
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	 );

public: // methods

	/// @brief Apply setup for membrane symmetry move
	virtual
	void
	apply( Pose & pose ); 

	/// @brief Get the name of this mover
	virtual std::string get_name() const; 

private: // helper methods

	/// @brief Helper Methods - add membrane residue (pre-symmdef)
	core::SSize
	setup_asymm_membrane( Pose & pose );
	
	/// @brief Register options from the command line
	void
	register_options();
	
	/// @brief Initialize option from the command line
	void
	init_from_cmd();

private: // data

	std::string symmdef_file_;

};


} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetry_hh
