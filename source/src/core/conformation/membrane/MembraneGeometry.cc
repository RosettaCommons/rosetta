// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/conformation/membrane/MembraneGeometry.cc
/// @brief    Base class for geometry of the membrane
///
/// @details  MembraneGeometry is an interface class for different membrane geometries.
///
/// @note     This object is a member of Conformation and should only be accessed using
///           pose.conformation().membrane_geometry().
///
/// @author   Hope Woods (hope.woods@vanderbilt.edu)

// Unit Headers
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Project Headers
#include <core/conformation/membrane/MembraneParams.hh>

// Package Headers
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>

static basic::Tracer TR( "core.conformation.membrane.MembraneGeometry" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace membrane {

// class Conformation;
// class MembraneInfo;


/// @brief Create MembraneInfo from initialized data
MembraneGeometry::MembraneGeometry(
	core::Real steepness
) :
	thickness_( 15 ),
	steepness_( steepness )
{}

/// @brief Create MembraneInfo from initialized data
MembraneGeometry::MembraneGeometry(
	core::Real steepness,
	core::Real thickness
) :
	thickness_( thickness ),
	steepness_( steepness )
{}

/// @brief Destructor
MembraneGeometry::~MembraneGeometry() {}


//thickness_vector returns a normalized vector that when the protein is tranformed into
//membrane coordinates, should be in the same direction as the x-axis
core::Vector
MembraneGeometry::thickness_vector( Conformation const & conf ) const {
	core::Size mem_rsd_num = conf.membrane_info()->membrane_rsd_num();
	core::Vector thickness = conf.residue( mem_rsd_num ).xyz( membrane::thickness );
	core::Vector new_x_axis = thickness - conf.membrane_info()->membrane_center( conf );
	new_x_axis.normalize(); //normalize here should not be a problem, valid coordinate frame checked in Conformation::check_valid_membrane
	return new_x_axis;
}

//normal_vector returns a normalized vector in the direction of the membrane normal
//If transformed into membrane coordinates this would be in the same direction as the z-axis
core::Vector
MembraneGeometry::normal_vector( Conformation const & conf ) const {
	core::Size mem_rsd_num = conf.membrane_info()->membrane_rsd_num();
	core::Vector normal = conf.residue( mem_rsd_num ).xyz( membrane::normal );
	core::Vector new_z_axis = normal - conf.membrane_info()->membrane_center( conf );
	new_z_axis.normalize(); //normalize here should not be a problem, valid coordinate frame checked in Conformation::check_valid_membrane
	return new_z_axis;
}

//binormal_vector returns a normalized vector that is normal to the normal_vector
//and the thickness_vector. In membrane coordinates this would be in the same direction
//as the y-axis
core::Vector
MembraneGeometry::binormal_vector( Conformation const & conf ) const {
	core::Vector thickness_v = thickness_vector( conf );
	core::Vector normal_v = normal_vector( conf );
	core::Vector new_y_axis = normal_v.cross( thickness_v );
	//normal_v and thickness_v should be orthogonal so no need to normalize
	return new_y_axis;
}

//returns corrected coordinate for the position of an atom with respect to
//membrane location and orientation
core::Real
MembraneGeometry::corrected_coordinate( core::Vector const & xyz, core::Vector const & axis ) const {
	core::Real coordinate = xyz.dot( axis );
	return coordinate;
}


//returns corrected xyz coordiantes of an atom with respect to membrane location and orientation
core::Vector
MembraneGeometry::corrected_xyz( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	Vector const & xyz( conf.residue( resnum ).atom( atomnum ).xyz() );
	core::Vector xyz_c = xyz - mem_cen;
	core::Vector x_axis = thickness_vector( conf );
	core::Vector y_axis = binormal_vector( conf );
	core::Vector z_axis = normal_vector( conf );
	core::Real x = corrected_coordinate( xyz_c, x_axis );
	core::Real y = corrected_coordinate( xyz_c, y_axis );
	core::Real z = corrected_coordinate( xyz_c, z_axis );
	core::Vector correct_xyz{ x, y, z };
	return correct_xyz;
}

//f1 and f2 are needed for calculating the derivative during minimization
//described in  H. Abe, W. Braun, T. Noguti, N. Go, Computers & Chemistry 1984 and also the last ~30 minutes of this youtube video from bootcamp: https://www.youtube.com/watch?v=j07ibj-fT1A.
//f1 returns ((r_alpha X atom_xyz)/|r_alpha - atom_xyz|)*dE/dr, where dE/dr is the derivative of the transition function with respect to the distance between atom_xyz and r_alpha.
//Here, r_alpha represents the point in space where the derivative of the transition functions depends on the distance between r_alpha and the atom's xyz position.
core::Vector
MembraneGeometry::f1( core::Vector atom_xyz, core::Vector r_alpha, core::Real deriv ) const {

	core::Vector d = r_alpha - atom_xyz;
	core::Real d_norm = d.length();
	if ( d_norm == Real(0.0) ) return { 0.0,0.0,0.0 };

	core::Real invd = 1.0 / d_norm; core::Vector f1 = r_alpha.cross( atom_xyz );
	return f1*invd*deriv;
}

//f2 returns ((r_alpha - atom_xyz)/|r_alpha - atom_xyz|)*dE/dr, where dE/dr is the derivative of the transition function with respect to the distance between atom_xyz and r_alpha.
core::Vector
MembraneGeometry::f2( core::Vector atom_xyz, core::Vector r_alpha, core::Real deriv ) const {

	core::Vector d = r_alpha - atom_xyz;
	core::Real d_norm = d.length();
	if ( d_norm == Real(0.0) ) return { 0.0,0.0,0.0 };

	core::Real invd = 1.0 / d_norm;
	return d * invd * deriv;
}

//f_imm1 is utilized in slab and bicelle geometry so we define once here
//n is steepness of hydrophobic -> hydrophillic transition (default = 10)
//transition function from Lazaridis. Effective Energy Function for Proteins in Lipid Membranes. Proteins. 2003

core::Real
MembraneGeometry::f_imm1( core::Real z_position ) const {

	core::Real n(0), pot(0), f(0), t(0);
	t = membrane_thickness();
	n = membrane_steepness();
	pot = std::abs( z_position )/t;
	f = std::pow( pot, n )/(1 + std::pow( pot, n ));
	return f;
}

core::Real
MembraneGeometry::f_imm1_deriv( core::Real z_position ) const {

	core::Real t = membrane_thickness();
	core::Real n = membrane_steepness();
	core::Real pot = std::abs( z_position )/t;
	core::Real numerator = n*pow(pot, n-1);
	core::Real denominator = t*(pow(1 + pow( pot, n ), 2));
	core::Real f_dz_position = numerator/denominator;
	return f_dz_position;
}


/// @brief Effective thickness of the membrane (default = 15)
core::Real
MembraneGeometry::membrane_thickness() const {
	return thickness_;
}

/// @brief Steepness of hydrophobic -> hydrophillic transition (default = 15)
core::Real
MembraneGeometry::membrane_steepness() const {
	return steepness_;
}


} // membrane
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::membrane::MembraneGeometry::save( Archive & arc ) const {
	arc( CEREAL_NVP( thickness_ ) ); // core::Real
	arc( CEREAL_NVP( steepness_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::membrane::MembraneGeometry::load( Archive & arc ) {
	arc( thickness_ ); // core::Real
	arc( steepness_ ); // core::Real
}
SAVE_AND_LOAD_SERIALIZABLE( core::conformation::membrane::MembraneGeometry );
CEREAL_REGISTER_TYPE( core::conformation::membrane::MembraneGeometry )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_membrane_MembraneGeometry )
#endif // SERIALIZATION

