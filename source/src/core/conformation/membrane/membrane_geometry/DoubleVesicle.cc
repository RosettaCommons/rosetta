// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/membrane_geometry/DoubleVesicle.cc
/// @brief Data describing the parameters of a vesicle
///
/// @details DoubleVesicle class contains the parameters of a vesicle and
///  the function to calculate the transition from the hydrophobic
///  environment of the bilayer to a hydrophilic environment.
///  In membrane coordinates the origin is located at the center of
///  the inner vesicle membrane.

/// @note This object is a member of Conformation and should only be accessed using
///            pose.conformation().membrane_geometry().
///
/// @author Hope Woods (hope.woods@vanderbilt.edu)

// Unit Headers
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/membrane_geometry/DoubleVesicle.hh>

// Project Headers
#include <core/conformation/membrane/MembraneParams.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>

// Package Headers
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

// C++ Headers
#include <string>

static basic::Tracer TR( "core.conformation.membrane.membrane_geometry.DoubleVesicle" );

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
namespace membrane_geometry {


DoubleVesicle::DoubleVesicle(
	core::Real steepness
) :
	MembraneGeometry( steepness )
{
	update_radii();
}

DoubleVesicle::DoubleVesicle(
	core::Real steepness,
	core::Real thickness
) :
	MembraneGeometry( steepness, thickness )
{
	update_radii();
}

DoubleVesicle::DoubleVesicle(
	core::Real steepness,
	core::Real thickness,
	core::Real inner_radius,
	core::Real distance
) :
	MembraneGeometry( steepness, thickness )
{
	set_inner_radius( inner_radius );
	set_distance( distance );
}

DoubleVesicle::DoubleVesicle(
	core::Real steepness,
	core::Real thickness,
	core::Real inner_radius,
	core::Real distance,
	AqueousPoreParametersOP aqueous_pore
) :
	MembraneGeometry( steepness, thickness, aqueous_pore )
{
	set_inner_radius( inner_radius );
	set_distance( distance );
}


/// @brief Destructor
DoubleVesicle::~DoubleVesicle() {}

MembraneGeometryOP DoubleVesicle::clone() const {
	return DoubleVesicleOP( new DoubleVesicle( *this ) );
}

/// @brief Generate a string representation of information represented by DoubleVesicle
void
DoubleVesicle::show() const {
	show( std::cout );
}

void
DoubleVesicle::show( std::ostream & output ) const {
	output << "MembraneGeometry: Information about the geometry of the membrane or membrane mimetic" << std::endl;
	output << "Outer Vesicle Radius: " << outer_radius_ << std::endl;
	output << "Inner Vesicle Radius: " << inner_radius_ << std::endl;
}

void
DoubleVesicle::update_radii( ) {
	if ( basic::options::option[ basic::options::OptionKeys::mp::geo::double_vesicle_inner_radius ].user() ) {
		set_inner_radius( basic::options::option[ basic::options::OptionKeys::mp::geo::double_vesicle_inner_radius ]() );
		TR << "setting inner vesicle radius from command line option:  " << inner_radius_ << std::endl;
	} else {
		TR << "setting vesicle inner_radius_ as default 1000" << std::endl;
		set_inner_radius( 1000 );
	}

	if ( basic::options::option[ basic::options::OptionKeys::mp::geo::double_vesicle_distance ].user() ) {
		set_distance( basic::options::option[ basic::options::OptionKeys::mp::geo::double_vesicle_distance ]());
	} else {
		set_distance(30);
	}

	set_outer_radius( inner_radius_ + ((membrane_thickness()*2)+distance_) );
	TR << "Setting outer radius: " << outer_radius_ << std::endl;
}

void
DoubleVesicle::set_distance( core::Real distance ) {
	if ( distance < 0 ) {
		TR.Fatal << "Cannot set distance between vesicle membranes as a negative number" << std::endl;
	} else if ( distance < membrane_thickness() ) {
		TR.Warning << "Distance between vesicle membranes is less than membrane thickness." << std::endl;
		TR.Warning << "May have unexpected behavior." << std::endl;
	}
	distance_ = distance;
	set_outer_radius( inner_radius_ + ((membrane_thickness()*2)+distance_) );
}

void
DoubleVesicle::set_outer_radius( core::Real outer_r ) {
	if ( outer_r < inner_radius_ ) {
		TR.Fatal << "Double Vesicle outer radius cannot be less than the inner radius." << std::endl;
	}
	outer_radius_ = outer_r;
}

void
DoubleVesicle::set_inner_radius( core::Real inner_r ) {
	if ( inner_r <= 0 ) {
		TR.Fatal << "Vesicle inner radius must be greater than 0" << std::endl;
	}
	inner_radius_ = inner_r;
}

core::Real
DoubleVesicle::get_outer_radius() const {
	return outer_radius_;
}

core::Real
DoubleVesicle::get_inner_radius() const {
	return inner_radius_;
}

core::Real
DoubleVesicle::get_distance() const {
	return distance_;
}

//DoubleVesicle transition function and helper functions
//xyz is the coordinates in space of the atom of interest
//n is steepness of hydrophobic -> hydrophillic transition (default = 15)
core::Real
DoubleVesicle::center( core::Vector xyz ) const {
	return pow( pow( xyz.x(), 2 ) + pow( xyz.y(), 2 ) + pow( xyz.z()+ inner_radius_ , 2 ), 0.5);
}

core::Real
DoubleVesicle::f_vesicle_membrane( Conformation const & conf, core::Vector xyz, core::Real radius ) const {
	core::Real c = center( xyz );
	core::Real j( f_thickness( conf, radius-c) );
	return j;
}

core::Real
DoubleVesicle::f_double_vesicle( Conformation const & conf, core::Vector xyz ) const {
	core::Real o = f_vesicle_membrane( conf, xyz, outer_radius_ );
	core::Real i = f_vesicle_membrane( conf, xyz, inner_radius_ );
	return o*i;
}


core::Real
DoubleVesicle::f_deriv( Conformation const & conf, core::Vector xyz, core::Real radius ) const {
	core::Real c = center( xyz );
	return f_thickness_deriv( conf, radius - c );
}



core::Real
DoubleVesicle::f_vesicle_deriv( Conformation const & conf, core::Vector xyz ) const {
	core::Real f_inner_deriv = f_deriv( conf, xyz, inner_radius_ );
	core::Real f_outer_deriv = f_deriv( conf, xyz, outer_radius_ );
	return f_inner_deriv*f_vesicle_membrane( conf, xyz, outer_radius_ ) + f_outer_deriv*f_vesicle_membrane( conf, xyz, inner_radius_ );
}
//returns the value of the transition function for membrane score functions
core::Real
DoubleVesicle::f_transition( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real f( f_double_vesicle( conf, xyz ) );

	return f_hydration( f, xyz );
}

core::Real
DoubleVesicle::f_transition_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	return f_vesicle_deriv( conf, xyz );
}

//returns electrostatic field due to the DoubleVesicle geometry
//xyz is the coordinates in space of the atom of interest
core::Real
DoubleVesicle::fa_elec_lipid( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	TR.Warning << "This term is only configured for Slab geometry. This term is set to 0.0 in this geometry" << std::endl;
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real f_elec_lipid( conf.membrane_info()->implicit_lipids()->f_elec_field( xyz.z() ) );
	return f_elec_lipid;

}

core::Real
DoubleVesicle::fa_elec_lipid_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {

	TR.Warning << "This term is only configured for Slab geometry. This term is set to 0.0 in this geometry" << std::endl;
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real f_elec_lipid_deriv( conf.membrane_info()->implicit_lipids()->f_elec_field_gradient( xyz.z() ) );
	return f_elec_lipid_deriv;

}

core::Vector
DoubleVesicle::r_alpha( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	const core::Vector normal = conf.membrane_info()->membrane_normal( conf );
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Vector vesicle_center = mem_cen - (inner_radius_*normal);
	core::Vector P = xyz - vesicle_center;
	core::Vector Q = (inner_radius_/P.length())*P;
	core::Vector R = Q + vesicle_center;
	return R;
}

core::Vector
DoubleVesicle::f_transition_f1( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real deriv(f_transition_deriv( conf, resnum, atomnum ) );
	core::Vector r(r_alpha( conf, resnum, atomnum ));

	core::Vector f1_dvesicle( f1( xyz, r, deriv));

	if ( !has_pore() ) {
		return f1_dvesicle;
	}

	core::Real f_dves( f_double_vesicle( conf, xyz ) );
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	core::Vector f1_p( f1_pore( f_dves, xyz, conf, resnum, atomnum) );
	return f1_dvesicle + f1_p;
}

core::Vector
DoubleVesicle::f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real deriv(f_transition_deriv( conf, resnum, atomnum ) );
	core::Vector r(r_alpha( conf, resnum, atomnum ));

	core::Vector f2_dvesicle( f2( xyz, r, deriv));

	if ( !has_pore() ) {
		return f2_dvesicle;
	}

	core::Real f_dves( f_double_vesicle( conf, xyz ) );
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	core::Vector f2_p( f2_pore( f_dves, xyz, conf, resnum, atomnum) );
	return f2_dvesicle + f2_p;
}

//returning string of name of geometry that was created
std::string
DoubleVesicle::geometry_string( ) const {
	std::string geometry_name = "double_vesicle";
	return geometry_name;
}

//return geometry enum
MP_GEOMETRY_TRANSITION
DoubleVesicle::geometry_enum() const {
	return MP_GEOMETRY_TRANSITION::DOUBLE_VESICLE;
}

} // geometry
} // membrane
} // conformation
} // core

