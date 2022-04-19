// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/membrane_geometry/Slab.hh
/// @brief Data describing the parameters of a flat membrane, the Imm1 representation
///
/// @details Slab class contains the parameters of the imm1 representation
//      and the function to calculate the transition from the hydrophobic
///      environment of the bicelle to a hydrophilic environment.
///
/// @note This object is a member of Conformation and should only be accessed using
///            pose.conformation().membrane_geometry().
///
/// @author Hope Woods (hope.woods@vanderbilt.edu)

#ifndef INCLUDED_core_conformation_membrane_geometry_Slab_hh
#define INCLUDED_core_conformation_membrane_geometry_Slab_hh

// Unit headers
#include <core/conformation/membrane/membrane_geometry/Slab.fwd.hh>

// Package Headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>

// Project Headers
#include <core/types.hh>

// Utility Headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace membrane {
namespace membrane_geometry {

class Slab : public MembraneGeometry {

public: // Constructors & Setup

	/// @brief
	/// @details
	Slab( core::Real steepness );

	/// @brief
	/// @details
	Slab( core::Real steepness, core::Real thickness );

	/// @brief Destructor
	~Slab();

	SlabOP
	clone() const;

	// @brief Generate a string representation of information represented by Slab
	void show() const override;

	void show( std::ostream & output ) const override;


public:
	//returns the value of the transition function for membrane score functions
	core::Real f_transition( Conformation const & conf, core::Size resnum, core::Size atomnum ) const override;

	core::Real f_transition_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	core::Vector r_alpha( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	core::Vector f_transition_f1( Conformation const & conf, core::Size resnum, core::Size atomnum ) const override;
	core::Vector f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const override;

	//returning string of name of geometry that was created
	std::string geometry_string() const override;

	//return geometry enum
	MP_GEOMETRY_TRANSITION geometry_enum() const override;

};
} // geometry
} // membrane
} // conformation
} // core

#endif //INCLUDED_core_conformation_membrane_geometry_Slab_hh
