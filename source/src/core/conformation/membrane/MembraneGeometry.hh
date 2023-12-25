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
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.fwd.hh>
#include <core/conformation/membrane/AqueousPoreParameters.fwd.hh>

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

	/// @brief
	/// @details
	MembraneGeometry(
		core::Real steepness,
		core::Real thickness,
		AqueousPoreParametersOP aqueous_pore
	);

	/// @brief Destructor
	~MembraneGeometry() override;

	virtual MembraneGeometryOP
	clone() const = 0;

	/// @brief Generate a string representation of information represented by this MembraneGeometry and send it to std::cout
	virtual void show() const = 0;

	/// @brief Generate a string representation of information represented by this MembraneGeometry
	virtual void show( std::ostream & output ) const = 0;


	/// @brief Does this protein have a water-filled pore?
	bool has_pore() const;

protected: //ensuring atom coordinates are correct with respect to membrane coordinates

	//thickness_vector returns a normalized vector that when the protein is tranformed into
	//membrane coordinates, should be in the same direction as the x-axis
	core::Vector thickness_vector( Conformation const & conf ) const;

	//normal_vector returns a normalized vector in the direction of the membrane normal (with respect to a flat membrane)
	//If transformed into membrane coordinates this would be in the same direction as the z-axis
	core::Vector normal_vector( Conformation const & conf ) const;

	//binormal_vector returns a normalized vector that is normal to the normal_vector
	//and the thickness_vector. In membrane coordinates this would be in the same direction
	//as the y-axis
	core::Vector binormal_vector( Conformation const & conf ) const;

	//returns corrected coordinate for the position of an atom with respect to
	//membrane location and orientation
	core::Real corrected_coordinate( core::Vector const & xyz, core::Vector const & x_axis ) const;



public:

	//returns corrected xyz coordiantes of an atom with respect to membrane location and orientation
	core::Vector corrected_xyz( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	// f_imm1 is utilized in slab and bicelle geometry so we define once here
	// transition function from Lazaridis. Effective Energy Function for Proteins in Lipid Membranes. Proteins. 2003
	core::Real f_imm1( core::Real z_position ) const;

	core::Real f_imm1_deriv( core::Real z_position ) const;

	// f_franklin is a transition function from:
	// Alfrod, Fleming, Fleming, Gray. Protein Structure Prediction and Design in a Biologically Realistic Implicit Membrane. 2020.
	// Tau and kappa are parameters that depend on the Lipid type, defined in ImplicitLipidInfo
	core::Real f_franklin( core::Real const z, core::Real tau, core::Real kappa ) const;

	core::Real f_franklin_gradient( core::Real const z, core::Real tau, core::Real kappa ) const;


	// returns transition value for flat membrane with no pore
	// if use_franklin == true returns f_franklin transition
	// else it returns f_imm1
	// This function does NOT take pore into account
	core::Real f_thickness( Conformation const & conf, core::Real const z ) const;

	//returns derivative of the transition value for a flat membrane with no pore
	core::Real f_thickness_deriv( Conformation const & conf, core::Real const z ) const;

	// Checks for existence of a pore
	// If there isn't a pore, returns the same transition function value as the input f_thk
	// If there is a pore, it calculates the pore transition function value (f_cavity) and
	// returns the composition of the pore and f_thk
	core::Real f_hydration( core::Real f_thk, numeric::xyzVector< core::Real > const & p ) const;

	//core::Real f_hydration_deriv_dz( numeric::xyzVector< core::Real > const & p, core::Real f_thk, core::Real f_thk_deriv_dz ) const;
	core::Real f_hydration_deriv_dz( numeric::xyzVector< core::Real > const & p, core::Real f_thk_deriv_dz ) const;

	//Pore transition function and derivative

	/// @brief Calculate the hydration of an atom based on its location relative to
	/// an aqueous pore or cavity
	core::Real
	f_cavity( numeric::xyzVector< core::Real > const & p ) const;

	/// @brief Calculate the derivative of f_cavity (without any r(x,y,z) dependence)
	core::Real
	f_cavity_gradient( core::Real const r ) const;

	// partial derivative with respect to x of the pore transition function
	core::Real f_cavity_dx( numeric::xyzVector< core::Real > const & p, core::Real f_thk ) const;

	// partial derivative with respect to y of the pore transition function
	core::Real f_cavity_dy( numeric::xyzVector< core::Real > const & p, core::Real f_thk ) const;

	// partial derivative with respect to z of the pore transition function
	core::Real f_cavity_dz( numeric::xyzVector< core::Real > const & p, core::Real f_thk ) const;

	/// @brief Calculate the location of an atom relative to the pore structure
	core::Real
	g_radius( numeric::xyzVector< core::Real > const & p ) const;

	// Derivative of g_radius with respect to z
	core::Real
	g_radius_gradient_dz( numeric::xyzVector< core::Real > const & p ) const;

	// Derivative of g_radius with respect to x
	core::Real
	g_radius_gradient_dx( numeric::xyzVector< core::Real > const & p ) const;

	// Derivative of g_radius with respect to y
	core::Real
	g_radius_gradient_dy( numeric::xyzVector< core::Real > const & p ) const;

	// r_alpha for calculating the derivative of the pore transition function with respect to x
	// See comment for f1 below for description of r_alpha
	core::Vector
	r_alpha_p_x( numeric::xyzVector< core::Real > const & xyz ) const;

	// r_alpha for calculating the derivative of the pore transition function with respect to y
	// See comment for f1 below for description of r_alpha
	core::Vector
	r_alpha_p_y( numeric::xyzVector< core::Real > const & xyz ) const;

	// r_alpha for calculating the derivative of the pore transition function with respect to z
	// See comment for f1 below for description of r_alpha
	core::Vector
	r_alpha_p_z( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;
	//r_alpha_p_z(  numeric::xyzVector< core::Real > const & xyz, const core::Vector mem_cen ) const;

	//Calculate the sum of the x and y portions of f1 vector for the pore
	core::Vector
	f1_pore( core::Real f_thk, numeric::xyzVector< core::Real > const & xyz, Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	//Calculate the sum of the x and y portions of f2 vector for the pore
	core::Vector
	f2_pore( core::Real f_thk, numeric::xyzVector< core::Real > const & xyz, Conformation const & conf, core::Size resnum, core::Size atomnum ) const;


public: // Information about the membrane geometry


	/// @brief Set membrane aqueous pore parameters
	void set_aqueous_pore_parameters( AqueousPoreParametersOP aqueous_pore );

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


	//f1 and f2 are needed for calculating the derivative during minimization
	//described in  H. Abe, W. Braun, T. Noguti, N. Go, Computers & Chemistry 1984 and also the last ~30 minutes of this youtube video from bootcamp: https://www.youtube.com/watch?v=j07ibj-fT1A.
	//f1 returns ((r_alpha X atom_xyz)/|r_alpha - atom_xyz|)*dE/dr, where dE/dr is the derivative of the transition function with respect to the distance between atom_xyz and r_alpha.
	//Here, r_alpha represents the point in space where the derivative of the transition functions depends on the distance between r_alpha and the atom's xyz position.
	core::Vector f1( core::Vector const & atom_xyz, core::Vector const & r_alpha, core::Real deriv ) const;

	//f2 returns ((r_alpha - atom_xyz)/|r_alpha - atom_xyz|)*dE/dr, where dE/dr is the derivative of the transition function with respect to the distance between atom_xyz and r_alpha.
	core::Vector f2( core::Vector const & atom_xyz, core::Vector const & r_alpha, core::Real deriv ) const;




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

	virtual core::Vector
	f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const = 0;

	//fa_elec_lipid calculates the charge times electrostatic potential due to the lipid layer
	virtual core::Real
	fa_elec_lipid( Conformation const & conf, core::Size resnum, core::Size atomnum ) const = 0;

	//fa_elec_lipid_deriv calculates the derivativ due to fa_elec_lipid
	virtual core::Real
	fa_elec_lipid_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const = 0;



private: // data

	// membrane thickness is half the bilayer in angstroms
	core::Real thickness_;

	// membrane steepness represents the steepness of the transition
	// from a hydrophobic envrionment to a hydrophilic environment.
	// Lazaridis. Effective Energy Function for Proteins in Lipid Membranes. Proteins. 2003
	core::Real steepness_;


	// aqueous pore parameters
	core::Real pore_transition_steepness_;
	core::conformation::membrane::AqueousPoreParametersOP pore_params_;


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



