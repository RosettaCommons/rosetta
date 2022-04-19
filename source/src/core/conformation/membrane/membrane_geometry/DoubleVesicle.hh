// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/membrane_geometry/DoubleVesicle.hh
/// @brief Data describing the parameters of a vesicle
///
/// @details DoubleVesicle class contains the parameters of a vesicle
///      and the function to calculate the transition from the hydrophobic
///      environment of the bilayer to a hydrophilic environment.
///      In membrane coordinates the origin is located at the center of
///      the outer vesicle membrane.
///
/// @note This object is a member of Conformation and should only be accessed using
///            pose.conformation().membrane_geometry().
///
/// @author Hope Woods (hope.woods@vanderbilt.edu)

#ifndef INCLUDED_core_conformation_membrane_geometry_DoubleVesicle_hh
#define INCLUDED_core_conformation_membrane_geometry_DoubleVesicle_hh

// Unit headers
#include <core/conformation/membrane/membrane_geometry/DoubleVesicle.fwd.hh>

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

class DoubleVesicle : public MembraneGeometry {

public: // Constructors & Setup

	/// @brief
	/// @details
	DoubleVesicle( core::Real steepness );

	/// @brief
	/// @details
	DoubleVesicle( core::Real steepness, core::Real thickness );

	DoubleVesicle( core::Real steepness, core::Real thickness, core::Real outer_radius, core::Real distance );

	/// @brief Destructor
	~DoubleVesicle();

	DoubleVesicleOP
	clone() const;

	// @brief Generate a string representation of information represented by Bicelle
	void show() const override;

	void show( std::ostream & output ) const override;

	void set_distance( core::Real distance );
	void set_outer_radius( core::Real outer_r );
	void set_inner_radius( core::Real inner_r );

	void update_radii();

	core::Real get_outer_radius() const;
	core::Real get_inner_radius() const;
	core::Real get_distance() const;

protected:
	core::Real center( core::Vector xyz ) const;
	core::Real f_vesicle_membrane( core::Vector xyz, core::Real radius ) const;
	core::Real f_double_vesicle( core::Vector xyz ) const;
	core::Real f_deriv( core::Vector xyz, core::Real radius ) const;
	core::Real f_vesicle_deriv( core::Vector xyz ) const;

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

private:

	core::Real outer_radius_; //distance from vesicle center to membrane center of outer vesicle membrane
	core::Real inner_radius_; //distance from vesicle center to membrane center of inner vesicle membrane
	core::Real distance_; //distance between outer edge of inner membrane and inner edge of outer membrane

};
} // geometry
} // membrane
} // conformation
} // core

#endif //INCLUDED_core_conformation_membrane_geometry_DoubleVesicle_hh
