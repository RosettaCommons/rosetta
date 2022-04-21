// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/membrane_geometry/Bicelle.hh
/// @brief Data describing the parameters of the bicelle
///
/// @details Bicelle class contains the parameters of the bicelle and
///  the function to calculate the transition from the hydrophobic
///  environment of the bicelle to a hydrophilic environment.
///
/// @note This object is a member of Conformation and should only be accessed using
///            pose.conformation().membrane_geometry().
///
/// @author Hope Woods (hope.woods@vanderbilt.edu)

#ifndef INCLUDED_core_conformation_membrane_geometry_Bicelle_hh
#define INCLUDED_core_conformation_membrane_geometry_Bicelle_hh

// Unit headers
#include <core/conformation/membrane/membrane_geometry/Bicelle.fwd.hh>

// Package Headers
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/Conformation.fwd.hh>

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


class Bicelle : public MembraneGeometry {


public: // Constructors & Setup

	/// @brief
	/// @details
	Bicelle( core::Real steepness, Conformation const & conf, core::Size membrane_pos );

	/// @brief
	/// @details
	Bicelle( core::Real steepness, core::Real thickness, Conformation const & conf, core::Size membrane_pos );

	/// @brief
	/// @details
	Bicelle( core::Real steepness, core::Real thickness, core::Real bicelle_inner_radius );

	/// @brief Destructor
	~Bicelle() override;

	BicelleOP
	clone() const;

	// @brief Generate a string representation of information represented by Bicelle
	void show() const override;

	void show( std::ostream & output ) const override;

public: // Information about bicelle, more detailed description of parameters found below where private member variables are declared

	void update_radii();

	//Sets protein_slice_diameter_, should be calculated by protein_slice_diameter_at_mem_cen
	//calls update_radii() to update bicelle_inner_radius and bicelle_outer_radius after protein_slice_diameter_ is set.
	void set_protein_slice_diameter( core::Real diameter );

	void set_inner_radius( core::Real inner_r);
	void set_outer_radius ( core::Real outer_r);
	void set_bicelle_edge_steepness( core::Real edge_steepness);

	core::Real protein_slice_diameter() const;

	core::Real bicelle_edge_steepness() const;

	core::Real bicelle_inner_radius() const;

	core::Real bicelle_outer_radius() const;

private: // default constructor

	Bicelle();

private:

	void update_edge_steepness();

	//This function calculates the distance between every CA atom within 3A of the membrane center plane.
	//In membrane coordinates the membrane center plane is the xy plane with the z axis as the membrane normal.
	//Instead of returning the absolute max distance found, it returns the 95th percentile of the CA distances.
	core::Real protein_slice_diameter_at_mem_cen( Conformation const & conf, core::Size membrane_pos ) const;


	void update_inner_radius();

	void update_outer_radius();

	// Bicelle transition function

	//transition function for bicelle edge
	core::Real h_bicelle( core::Vector xyz, const core::Vector mem_cen ) const;

	//combine h_bicelle and f_imm1
	core::Real f_bicelle( core::Vector xyz, core::Real z_depth, const core::Vector mem_cen ) const;

	//derivative of h_bicelle with respect to r
	core::Real h_bicelle_deriv_wrt_r( core::Vector xyz, const core::Vector mem_cen ) const;

	//derivative of f_bicelle
	core::Real f_bicelle_deriv( core::Vector xyz, core::Real z_depth, const core::Vector mem_cen ) const;

public:

	core::Real f_transition( Conformation const & conf, core::Size resnum, core::Size atomnum ) const override;

	core::Real f_transition_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	core::Vector r_alpha( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;
	core::Vector r_alpha_m( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	core::Vector f_transition_f1( Conformation const & conf, core::Size resnum, core::Size atomnum ) const override;
	core::Vector f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const override;


	core::Real f_transition_deriv_m( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;




	//returning string of name of geometry that was created
	std::string geometry_string() const override;

	//return geometry enum
	MP_GEOMETRY_TRANSITION geometry_enum() const override;

private:

	// /*                 ***************************************
	//                  **                                      **
	//                 **                                        **
	//                **                       outer_radius       **
	//               **                     *----------------------**
	//              **                      *-------------------   **
	//               **                      inner_radius      |  **
	//                **                                       | **
	//                 **                                      |**
	//                   ***************************************
	// bicelle_inner_radius is either defined by the user or set automatically
	// the bicelle_outer_radius is defined as the length of the bicelle_inner_radius + membrane_thickness.
	// membrane_thickness is actally the thickness of one leaflet (half the membrane)

	// The bicelle_edge_steepness is the steepness of the transition from a hydrophbic
	// environment to a hydrophilic environment on the curved edge of the bicelle;
	// where the double astricks (**) are in the above drawing.

	core::Real bicelle_inner_radius_ = 0.0;
	core::Real bicelle_outer_radius_ = 0.0;

	core::Real bicelle_edge_steepness_ = 0.0;

	//longest distance between CA atoms from a protein slice in the membrane center
	core::Real protein_slice_diameter_ = 0.0;

};
} // geometry
} // membrane
} // conformation
} // core

#endif //INCLUDED_core_conformation_membrane_geometry_Bicelle_hh
