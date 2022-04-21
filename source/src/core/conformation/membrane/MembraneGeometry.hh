// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/conformation/membrane/MembraneGeometry.hh
/// @brief    Base class for geometry of the membrane
///
/// @details  MembraneGeometry is a container object that describes the geometry of the membrane
///
/// @note     This object is a member of Conformation and should only be accessed using
///           pose.conformation().membrane_geometry().
///
/// @author   Hope Woods (hope.woods@vanderbilt.edu)

#ifndef INCLUDED_core_conformation_membrane_MembraneGeometry_hh
#define INCLUDED_core_conformation_membrane_MembraneGeometry_hh

// Unit headers
#include <core/conformation/membrane/MembraneGeometry.fwd.hh>

// Package Headers
#include <core/conformation/Conformation.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/VirtualBase.hh>

// C++ Headers
#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace membrane {


/// @brief Data describing the geometry of the membrane
class MembraneGeometry : public utility::VirtualBase {

public: // Constructors & Setup

	/// @brief Create a default version of MembraneGeometry (DONT USE)
	MembraneGeometry() = delete;


	/// @brief
	/// @details
	MembraneGeometry(
		core::Real steepness
	);

	/// @brief
	/// @details
	MembraneGeometry(
		core::Real steepness,
		core::Real thickness
	);

	/// @brief Destructor
	~MembraneGeometry() override;

	/// @brief Generate a string representation of information represented by this MembraneGeometry and send it to std::cout
	virtual void show() const = 0;

	/// @brief Generate a string representation of information represented by this MembraneGeometry
	virtual void show( std::ostream & output ) const = 0;

protected:

	//f_imm1 is utilized in slab and bicelle geometry so we define once here
	core::Real f_imm1( core::Real z_position ) const;


	core::Real f_imm1_deriv( core::Real z_position ) const;

	//thickness_vector returns a normalized vector that when the protein is tranformed into
	//membrane coordinates, should be in the same direction as the x-axis
	core::Vector thickness_vector( Conformation const & conf ) const;

	//binormal_vector returns a normalized vector that is normal to the normal_vector
	//and the thickness_vector. In membrane coordinates this would be in the same direction
	//as the y-axis
	core::Vector binormal_vector( Conformation const & conf ) const;

	//returns corrected coordinate for the position of an atom with respect to
	//membrane location and orientation
	core::Real corrected_coordinate( core::Vector const & xyz, core::Vector const & x_axis ) const;



public: // Information about the membrane geometry

	//normal_vector returns a normalized vector in the direction of the membrane normal (with respect to a flat membrane)
	//If transformed into membrane coordinates this would be in the same direction as the z-axis
	core::Vector normal_vector( Conformation const & conf ) const;


	//returning string of name of geometry that was created
	virtual std::string
	geometry_string() const = 0;

	//return geometry enum
	virtual MP_GEOMETRY_TRANSITION
	geometry_enum() const =0;

	/// @brief Effective thickness of the membrane (default = 15)
	virtual core::Real membrane_thickness() const;

	/// @brief Steepness of hydrophobic -> hydrophillic transition (defualt = 10)
	virtual core::Real membrane_steepness() const;


	//returns corrected xyz coordiantes of an atom with respect to membrane location and orientation
	core::Vector corrected_xyz( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	//f1 and f2 are needed for calculating the derivative during minimization
	//described in  H. Abe, W. Braun, T. Noguti, N. Go, Computers & Chemistry 1984 and also the last ~30 minutes of this youtube video from bootcamp: https://www.youtube.com/watch?v=j07ibj-fT1A.
	//f1 returns ((r_alpha X atom_xyz)/|r_alpha - atom_xyz|)*dE/dr, where dE/dr is the derivative of the transition function with respect to the distance between atom_xyz and r_alpha.
	//Here, r_alpha represents the point in space where the derivative of the transition functions depends on the distance between r_alpha and the atom's xyz position.
	core::Vector f1( core::Vector atom_xyz, core::Vector r_alpha, core::Real deriv ) const;

	//f2 returns ((r_alpha - atom_xyz)/|r_alpha - atom_xyz|)*dE/dr, where dE/dr is the derivative of the transition function with respect to the distance between atom_xyz and r_alpha.
	core::Vector f2( core::Vector atom_xyz, core::Vector r_alpha, core::Real deriv ) const;


	// The following functions are pure virtual functions,
	// derived classes are responsible for calculating f_transition, f_transition_f1, and f_transition_f2.
	// These are all called in energy methods FaMPSolvEnergy and FaMPEnvEnergy.

	//f_transition maps the transition of the hydrophobic membrane environment to an aqueous environment.
	virtual core::Real
	f_transition( Conformation const & conf, core::Size resnum, core::Size atomnum ) const = 0;

	//Functions to calculate f1 and f2 for derivatives of transition functions. See the description above for f1 and f2.
	//For some transition function geometries, we have to calculate multiple f1 and f2 values for all partial derivatives.
	//f_transition_f1 should return the sum of f1 for all partial derivatives
	virtual core::Vector
	f_transition_f1( Conformation const & conf, core::Size resnum, core::Size atomnum ) const = 0;

	//f_transition_f1 should return the sum of f1 for all partial derivatives
	virtual core::Vector
	f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const = 0;

private: // data

	// membrane thickness is half the bilayer in angstroms
	core::Real thickness_;

	// membrane steepness represents the steepness of the transition
	// from a hydrophobic envrionment to a hydrophilic environment.
	// Lazaridis. Effective Energy Function for Proteins in Lipid Membranes. Proteins. 2003
	core::Real steepness_;





#ifdef    SERIALIZATION
public:
	friend class cereal::access;
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // membrane
} // conformation
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_membrane_MembraneGeometry )
#endif // SERIALIZATION


#endif // INCLUDED_core_conformation_membrane_MembraneGeometry_hh



