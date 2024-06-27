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
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/membrane/AqueousPoreParameters.hh>

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

/// @brief Create MembraneInfo from initialized data
MembraneGeometry::MembraneGeometry(
	core::Real steepness
) :
	thickness_( 15 ),
	steepness_( steepness ),
	pore_transition_steepness_( 10.0 ),
	pore_params_()
{}

/// @brief Create MembraneInfo from initialized data
MembraneGeometry::MembraneGeometry(
	core::Real steepness,
	core::Real thickness
) :
	thickness_( thickness ),
	steepness_( steepness ),
	pore_transition_steepness_( 10.0 ),
	pore_params_()
{}

/// @brief Create MembraneGeometry from initialized data
MembraneGeometry::MembraneGeometry(
	core::Real steepness,
	core::Real thickness,
	AqueousPoreParametersOP aqueous_pore
) :
	thickness_( thickness ),
	steepness_( steepness ),
	pore_transition_steepness_( 10.0 ),
	pore_params_()
{
	set_aqueous_pore_parameters( aqueous_pore );
}


/// @brief Destructor
MembraneGeometry::~MembraneGeometry() {}


/// @brief Are we accommodating the aqueous pore?
bool
MembraneGeometry::has_pore() const {
	if ( pore_params_ == nullptr ) {
		return false;
	} else {
		return true;
	}
}

/// @brief Set membrane aqueous pore parameters
void
MembraneGeometry::set_aqueous_pore_parameters(
	AqueousPoreParametersOP aqueous_pore
) {
	pore_params_ = aqueous_pore;
}

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
MembraneGeometry::f1( core::Vector const & atom_xyz, core::Vector const & r_alpha, core::Real deriv ) const {

	core::Vector d = r_alpha - atom_xyz;
	core::Real d_norm = d.length();
	if ( d_norm == Real(0.0) ) return { 0.0,0.0,0.0 };

	core::Real invd = 1.0 / d_norm;
	core::Vector f1 = r_alpha.cross( atom_xyz );
	return f1*invd*deriv;
}

//f2 returns ((r_alpha - atom_xyz)/|r_alpha - atom_xyz|)*dE/dr, where dE/dr is the derivative of the transition function with respect to the distance between atom_xyz and r_alpha.
core::Vector
MembraneGeometry::f2( core::Vector const & atom_xyz, core::Vector const & r_alpha, core::Real deriv ) const {

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
	core::Real pot_n = std::pow( pot, n );

	f = pot_n/(1 + pot_n);
	return f;
}

core::Real
MembraneGeometry::f_imm1_deriv( core::Real z_position ) const {

	core::Real t = membrane_thickness();
	core::Real n = membrane_steepness();
	core::Real pot = std::abs( z_position )/t;
	core::Real numerator = n*std::pow(pot, n-1);
	core::Real pot_n_1 = 1 + std::pow( pot, n );
	core::Real denominator = t*(pot_n_1*pot_n_1);
	debug_assert( denominator != 0 );
	core::Real f_dz_position = numerator/denominator;
	return f_dz_position;
}


//f_franklin is a transition function from:
//Alfrod, Fleming, Fleming, Gray. Protein Structure Prediction and Design in a Biologically Realistic Implicit Membrane. 2020.
//Tau and kappa are parameters that depend on the Lipid type, defined in ImplicitLipidInfo
core::Real
MembraneGeometry::f_franklin( core::Real const z, core::Real tau, core::Real kappa ) const {
	core::Real abs_z( std::abs(z) );
	core::Real d( 1 + (tau*std::exp(-kappa*abs_z)));
	debug_assert( d != 0 );
	return 1/d;
}

core::Real
MembraneGeometry::f_franklin_gradient( core::Real const z, core::Real tau, core::Real kappa ) const {
	core::Real abs_z( std::abs(z) );
	core::Real exp_bz( std::exp(-kappa*abs_z) );
	core::Real denom = 1 + tau*exp_bz;
	debug_assert( denom != 0 );
	core::Real quotient( ( tau*kappa*exp_bz )/( denom * denom ) );
	return quotient;
}

// Checks for existence of a pore
// If there isn't a pore, returns the same transition function value as the input f_thk
// If there is a pore, it calculates the pore transition function value (f_cavity) and
// returns the composition of the pore and f_thk
core::Real
MembraneGeometry::f_hydration( core::Real f_thk, numeric::xyzVector< core::Real > const & p ) const {
	if ( !has_pore() ) {
		return f_thk;
	} else {
		core::Real f_cav( f_cavity( p ));
		return f_thk + f_cav - (f_thk*f_cav);
	}

}

core::Real
//MembraneGeometry::f_hydration_deriv_dz( numeric::xyzVector< core::Real > const & p, core::Real f_thk, core::Real f_thk_deriv_dz ) const {
MembraneGeometry::f_hydration_deriv_dz( numeric::xyzVector< core::Real > const & p, core::Real f_thk_deriv_dz ) const {
	core::Real f_cav( f_cavity( p ));

	return (f_thk_deriv_dz*(1-f_cav));
}


//returns transition value for flat membrane
//if use_franklin == true returns f_franklin transition
//else it returns f_imm1
//This function does NOT take pore into account
core::Real
MembraneGeometry::f_thickness( Conformation const & conf, core::Real const z ) const {

	core::Real f_thk;
	if ( conf.membrane_info()->use_franklin() ) {
		core::Real tau( conf.membrane_info()->implicit_lipids()->water_pseudo_thickness() );
		core::Real kappa( conf.membrane_info()->implicit_lipids()->water_steepness() );
		f_thk = f_franklin( z, tau, kappa) ;
	} else {
		f_thk = f_imm1( z );
	}

	return f_thk;
}


//returns derivative of the transition value for a flat membrane with no pore
core::Real
MembraneGeometry::f_thickness_deriv( Conformation const & conf, core::Real const z ) const {

	core::Real f_thk_deriv_dz;
	if ( conf.membrane_info()->use_franklin() ) {
		core::Real tau( conf.membrane_info()->implicit_lipids()->water_pseudo_thickness() );
		core::Real kappa( conf.membrane_info()->implicit_lipids()->water_steepness() );
		f_thk_deriv_dz = f_franklin_gradient( z, tau, kappa);
	} else {
		f_thk_deriv_dz = f_imm1_deriv( z );
	}

	return f_thk_deriv_dz;
}


// partial derivative with respect to x of the pore transition function
core::Real
MembraneGeometry::f_cavity_dx( numeric::xyzVector< core::Real > const & p, core::Real f_thk  ) const {
	core::Real r( g_radius( p ));
	core::Real dfcav_dg( f_cavity_gradient( r ) );

	core::Real dg_dx( g_radius_gradient_dx( p ));

	return dfcav_dg * dg_dx * ( 1 - f_thk );
}

core::Real
MembraneGeometry::f_cavity_dy( numeric::xyzVector< core::Real > const & p, core::Real f_thk  ) const {
	core::Real r( g_radius( p ));
	core::Real dfcav_dg( f_cavity_gradient( r ) );
	core::Real dg_dy( g_radius_gradient_dy( p ));

	return dfcav_dg * dg_dy * ( 1 - f_thk );
}

core::Real
MembraneGeometry::f_cavity_dz( numeric::xyzVector< core::Real > const & p, core::Real f_thk  ) const {
	core::Real r( g_radius( p ));
	core::Real dfcav_dg( f_cavity_gradient( r ) );
	core::Real dg_dz( g_radius_gradient_dz( p ));

	return dfcav_dg * dg_dz * ( 1 - f_thk );
}

/// @brief Calculate the hydration of an atom based on its location relative to
/// an aqueous pore or cavity
core::Real
MembraneGeometry::f_cavity( numeric::xyzVector< core::Real > const & p ) const {
	core::Real radius( g_radius(p) );
	core::Real r_n( std::pow( radius, pore_transition_steepness_ ) );
	core::Real quotient( r_n / (1+r_n) );
	return 1-quotient;
}

/// @brief Calculate the derivative of f_cavity (without any r(x,y,z) dependence)
core::Real
MembraneGeometry::f_cavity_gradient( core::Real const r ) const {
	core::Real top( pore_transition_steepness_*std::pow( r, pore_transition_steepness_-1 ) );
	core::Real r_trans_1 = std::pow( r, pore_transition_steepness_) + 1;
	core::Real bottom( r_trans_1*r_trans_1 );
	debug_assert( bottom != 0 );
	return -top/bottom;
}

/// @brief Calculate the location of an atom relative to the pore structure
core::Real
MembraneGeometry::g_radius( numeric::xyzVector< core::Real > const & p ) const {

	core::Real pore_center_x( pore_params_->pore_center_x( p.z() ) );
	core::Real pore_center_y( pore_params_->pore_center_y( p.z() ) );
	core::Real pore_minor_radius( pore_params_->pore_minor_radius( p.z() ) );
	core::Real pore_major_radius( pore_params_->pore_major_radius( p.z() ) );
	numeric::MathMatrix< core::Real > rotation( pore_params_->pore_rotation( p.z() ) );

	debug_assert( pore_minor_radius != 0 );
	debug_assert( pore_major_radius != 0 );

	core::Real abs_x_dis = std::abs(p.x() - pore_center_x);
	core::Real abs_y_dis = std::abs(p.y() - pore_center_y);

	core::Real lhs = (abs_x_dis*rotation(0,0) - abs_y_dis*rotation(0,1))/pore_major_radius;
	core::Real rhs = (abs_x_dis*rotation(1,0) - abs_y_dis*rotation(0,0))/pore_minor_radius;

	core::Real result = lhs*lhs + rhs*rhs;

	return result;
}

// Derivative of g_radius with respect to z
core::Real
MembraneGeometry::g_radius_gradient_dz( numeric::xyzVector< core::Real > const & p ) const {

	numeric::MathMatrix< core::Real > rotation( pore_params_->pore_rotation( p.z() ) );
	core::Real c( rotation(0,0) ); //cos_theta
	core::Real s( rotation(1,0) ); //sin_theta

	core::Real xo( pore_params_->pore_center_x( p.z() ) );
	core::Real yo( pore_params_->pore_center_y( p.z() ) );
	core::Real a( pore_params_->pore_major_radius( p.z() ) );
	core::Real b( pore_params_->pore_minor_radius( p.z() ) );
	debug_assert( a != 0 );
	debug_assert( b != 0 );

	core::Real d_theta( pore_params_->pore_rotation_deriv( p.z() ) );
	core::Real d_xo( pore_params_->pore_center_x_deriv( p.z() ) );
	core::Real d_yo( pore_params_->pore_center_y_deriv( p.z() ) );
	core::Real d_a( pore_params_->pore_major_radius_deriv( p.z() ) );
	core::Real d_b( pore_params_->pore_minor_radius_deriv( p.z() ) );

	core::Real x_dis = p.x()-xo;
	core::Real y_dis = p.y()-yo;
	core::Real abs_x_dis = std::abs(x_dis);
	core::Real abs_y_dis = std::abs(y_dis);
	core::Real sign_x_dis = std::copysign(1.0, x_dis);
	core::Real sign_y_dis = std::copysign(1.0, y_dis);

	core::Real h = (abs_x_dis*c + abs_y_dis*s);
	core::Real j = (abs_x_dis*s - abs_y_dis*c);

	core::Real dhdz = (sign_x_dis*(-d_xo * c) + abs_x_dis*(-s)*d_theta) + ( sign_y_dis*(-d_yo * s) + abs_y_dis*(+c)*d_theta);
	core::Real djdz = (sign_x_dis*(-d_xo * s) + abs_x_dis*(+c)*d_theta) - ( sign_y_dis*(-d_yo * c) + abs_y_dis*(-s)*d_theta);

	core::Real dlhs_dz = 2*((h/a)*((dhdz/a) - h*d_a/(a*a) ));
	core::Real drhs_dz = 2*((j/b)*((djdz/b) - j*d_b/(b*b) ));

	core::Real dgdz = dlhs_dz + drhs_dz;

	return dgdz;
}

// Derivative of g_radius with respect to x
core::Real
MembraneGeometry::g_radius_gradient_dx( numeric::xyzVector< core::Real > const & p ) const {

	core::Real c( pore_params_->pore_rotation( p.z() )(0,0) ); //cos_theta
	core::Real s( pore_params_->pore_rotation( p.z() )(1,0) ); //sin_theta
	core::Real ns( pore_params_->pore_rotation( p.z() )(0,1) ); //-sin_theta
	core::Real xo( pore_params_->pore_center_x( p.z() ) );
	core::Real yo( pore_params_->pore_center_y( p.z() ) );
	core::Real a( pore_params_->pore_major_radius( p.z() ) );
	core::Real b( pore_params_->pore_minor_radius( p.z() ) );
	debug_assert( a != 0 );
	debug_assert( b != 0 );

	core::Real abs_x_dis = std::abs(p.x() - xo);
	core::Real abs_y_dis = std::abs(p.y() - yo);

	core::Real dgdx( 2 * ( c/(a*a) * ( c*abs_x_dis - ns*abs_y_dis ) + s/(b*b) * ( s*abs_x_dis - c*abs_y_dis) ) );

	return dgdx;
}

// Derivative of g_radius with respect to y
core::Real
MembraneGeometry::g_radius_gradient_dy( numeric::xyzVector< core::Real > const & p ) const {

	core::Real c( pore_params_->pore_rotation( p.z() )(0,0) ); //cos_theta
	core::Real s( pore_params_->pore_rotation( p.z() )(1,0) ); //sin_theta
	core::Real ns( pore_params_->pore_rotation( p.z() )(0,1) ); //-sin_theta
	core::Real xo( pore_params_->pore_center_x( p.z() ) );
	core::Real yo( pore_params_->pore_center_y( p.z() ) );
	core::Real a( pore_params_->pore_major_radius( p.z() ) );
	core::Real b( pore_params_->pore_minor_radius( p.z() ) );
	debug_assert( a != 0 );
	debug_assert( b != 0 );

	core::Real abs_x_dis = std::abs(p.x() - xo);
	core::Real abs_y_dis = std::abs(p.y() - yo);

	core::Real dgdy( 2 * ( -ns/(a*a) * ( c*abs_x_dis - ns*abs_y_dis ) + -c/(b*b) * ( s*abs_x_dis - c*abs_y_dis) ) );

	return dgdy;
}

// r_alpha for calculating the derivative of the pore transition function with respect to x
// See comment for f1 below for description of r_alpha
core::Vector
MembraneGeometry::r_alpha_p_x( numeric::xyzVector< core::Real > const & xyz ) const {
	if ( ! has_pore() ) {
		core::Vector center = {0.0, 0.0, 0.0};
		return center;
	}
	core::Real pore_center_x( pore_params_->pore_center_x( xyz.z() ) );
	core::Vector center = { pore_center_x, xyz.y(), xyz.z() };
	return center;
}


// r_alpha for calculating the derivative of the pore transition function with respect to y
// See comment for f1 below for description of r_alpha
core::Vector
MembraneGeometry::r_alpha_p_y( numeric::xyzVector< core::Real > const & xyz ) const {
	if ( !has_pore() ) {
		core::Vector center = {0.0, 0.0, 0.0};
		return center;
	}
	core::Real pore_center_y( pore_params_->pore_center_y( xyz.z() ) );
	core::Vector center = { xyz.x(), pore_center_y, xyz.z() };
	return center;
}


// r_alpha for calculating the derivative of the pore transition function with respect to z
// See comment for f1 below for description of r_alpha
core::Vector
MembraneGeometry::r_alpha_p_z( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Real z_depth = conf.membrane_info()->atom_z_position( conf, resnum, atomnum );
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	const core::Vector normal = conf.membrane_info()->membrane_normal( conf );
	core::Vector const & xyz( conf.residue( resnum ).atom( atomnum ).xyz() );
	core::Vector proj_i = mem_cen + z_depth * normal;
	core::Vector i_ip = proj_i - xyz;
	return ( mem_cen - i_ip );
}
// MembraneGeometry::r_alpha_p_z( numeric::xyzVector< core::Real > const & xyz, const core::Vector mem_cen ) const {
//  core::Vector r_alpha_z = { mem_cen.x(), mem_cen.y(), xyz.z() };
//  return r_alpha_z;
// }


//Calculate the sum of the x and y portions of f1 vector for the pore
core::Vector
MembraneGeometry::f1_pore( core::Real f_thk, numeric::xyzVector< core::Real > const & xyz, Conformation const & conf, core::Size resnum, core::Size atomnum ) const {

	core::Real deriv_x( f_cavity_dx(xyz, f_thk) );
	core::Vector r_x( r_alpha_p_x( xyz ) );
	core::Vector f1_x( f1( xyz, r_x, deriv_x) );

	core::Real deriv_y( f_cavity_dy(xyz, f_thk) );
	core::Vector r_y( r_alpha_p_y( xyz ) );
	core::Vector f1_y( f1( xyz, r_y, deriv_y) );

	core::Real deriv_z( f_cavity_dz(xyz, f_thk) );
	core::Vector r_z( r_alpha_p_z( conf, resnum, atomnum ) );
	core::Vector f1_z( f1( xyz, r_z, deriv_z) );


	return f1_x + f1_y + f1_z;
}

//Calculate the sum of the x and y portions of f2 vector for the pore
core::Vector
MembraneGeometry::f2_pore( core::Real f_thk, numeric::xyzVector< core::Real > const & xyz, Conformation const & conf, core::Size resnum, core::Size atomnum ) const {

	core::Real deriv_x( f_cavity_dx(xyz, f_thk) );
	core::Vector r_x( r_alpha_p_x( xyz ) );
	core::Vector f2_x( f2( xyz, r_x, deriv_x) );

	core::Real deriv_y( f_cavity_dy(xyz, f_thk) );
	core::Vector r_y( r_alpha_p_y( xyz ) );
	core::Vector f2_y( f2( xyz, r_y, deriv_y) );

	core::Real deriv_z( f_cavity_dz(xyz, f_thk) );
	core::Vector r_z( r_alpha_p_z( conf, resnum, atomnum ) );
	core::Vector f2_z( f2( xyz, r_z, deriv_z) );


	return f2_x + f2_y + f2_z;
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
	arc( CEREAL_NVP( pore_transition_steepness_ ) ); // core::Real
	arc( CEREAL_NVP( pore_params_ ) ); // core::conformation::membrane::AqueousPoreParametersOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::membrane::MembraneGeometry::load( Archive & arc ) {
	arc( thickness_ ); // core::Real
	arc( steepness_ ); // core::Real
	arc( pore_transition_steepness_ ); // core::Real
	arc( pore_params_ ); // core::conformation::membrane::AqueousPoreParametersOP
}
SAVE_AND_LOAD_SERIALIZABLE( core::conformation::membrane::MembraneGeometry );
CEREAL_REGISTER_TYPE( core::conformation::membrane::MembraneGeometry )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_membrane_MembraneGeometry )
#endif // SERIALIZATION

