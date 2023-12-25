// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/membrane_geometry/Vesicle.cc
/// @brief Data describing the parameters of a vesicle
///
/// @details Vesicle class contains the parameters of a vesicle and
///  the function to calculate the transition from the hydrophobic
///  environment of the bilayer to a hydrophilic environment.
///
/// @note This object is a member of Conformation and should only be accessed using
///            pose.conformation().membrane_geometry().
///
/// @author Hope Woods (hope.woods@vanderbilt.edu)

// Unit Headers
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/membrane_geometry/Vesicle.hh>

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

static basic::Tracer TR( "core.conformation.membrane.membrane_geometry.Vesicle" );

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


Vesicle::Vesicle(
	core::Real steepness
) :
	MembraneGeometry( steepness )
{
	update_radius();
}

Vesicle::Vesicle(
	core::Real steepness,
	core::Real thickness
) :
	MembraneGeometry( steepness, thickness )
{
	update_radius();
}

Vesicle::Vesicle(
	core::Real steepness,
	core::Real thickness,
	core::Real radius
) :
	MembraneGeometry( steepness, thickness ),
	radius_( radius )
{}

Vesicle::Vesicle(
	core::Real steepness,
	core::Real thickness,
	core::Real radius,
	AqueousPoreParametersOP aqueous_pore
) :
	MembraneGeometry( steepness, thickness, aqueous_pore),
	radius_( radius )
{}

/// @brief Destructor
Vesicle::~Vesicle() {}

MembraneGeometryOP Vesicle::clone() const {
	return VesicleOP( new Vesicle( *this ) );
}

/// @brief Generate a string representation of information represented by Vesicle
void
Vesicle::show() const {
	show( std::cout );
}

void
Vesicle::show( std::ostream & output ) const {
	output << "MembraneGeometry: Information about the geometry of the membrane or membrane mimetic" << std::endl;
	output << "Vesicle Radius: " << radius_ << std::endl;
}

void
Vesicle::update_radius( ) {
	if ( basic::options::option[ basic::options::OptionKeys::mp::geo::vesicle_radius ].user() ) {
		set_radius(basic::options::option[ basic::options::OptionKeys::mp::geo::vesicle_radius ]());
		TR << "Setting vesicle radius from command line option:  " << radius_ << std::endl;
	} else {
		TR << "Setting vesicle radius_ as default 100" << std::endl;
		set_radius( 100 );
	}
}

void
Vesicle::set_radius( core::Real radius ) {
	radius_ = radius;
}

core::Real
Vesicle::get_radius() const {
	return radius_;
}

//Vesicle transition function and helper functions
//xyz is the coordinates in space of the atom of interest
//n is steepness of hydrophobic -> hydrophillic transition (default = 15)
core::Real
Vesicle::f_vesicle( Conformation const & conf, core::Vector xyz ) const {

	core::Real center = pow( pow( xyz.x(), 2 ) + pow( xyz.y(), 2 ) + pow( xyz.z()+ radius_ , 2 ), 0.5);
	core::Real f( f_thickness( conf, std::abs(radius_) - center) );

	return f;
}

core::Real
Vesicle::f_vesicle_deriv( Conformation const & conf, core:: Vector xyz) const {

	core::Real center = pow( pow( xyz.x(), 2 ) + pow( xyz.y(), 2 ) + pow( xyz.z()+ radius_ , 2 ), 0.5);
	core::Real df_dc( f_thickness_deriv( conf, std::abs(radius_) - center) );

	return df_dc;
}

//returns the value of the transition function for membrane score functions
core::Real
Vesicle::f_transition( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real f( f_vesicle( conf, xyz ) );

	return f_hydration( f, xyz );
}

core::Real
Vesicle::f_transition_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real f_ves_deriv( f_vesicle_deriv( conf, xyz) );
	if ( !has_pore() ) {
		return f_ves_deriv;
	} else {
		return f_hydration_deriv_dz( xyz, f_ves_deriv );
	}
}

//returns electrostatic field due to the Vesicle geometry
//xyz is the coordinates in space of the atom of interest
core::Real
Vesicle::fa_elec_lipid( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	TR.Warning << "This term is only configured for Slab geometry. This term is set to 0.0 in this geometry" << std::endl;
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real f_elec_lipid( conf.membrane_info()->implicit_lipids()->f_elec_field( xyz.z() ) );
	return f_elec_lipid;

}

core::Real
Vesicle::fa_elec_lipid_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {

	TR.Warning << "This term is only configured for Slab geometry. This term is set to 0.0 in this geometry" << std::endl;
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real f_elec_lipid_deriv( conf.membrane_info()->implicit_lipids()->f_elec_field_gradient( xyz.z() ) );
	return f_elec_lipid_deriv;

}

core::Vector
Vesicle::r_alpha( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	const core::Vector normal = conf.membrane_info()->membrane_normal( conf );
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Vector vesicle_center = mem_cen - (radius_*normal);
	core::Vector P = xyz - vesicle_center;
	core::Vector Q = (radius_/P.length())*P;
	core::Vector R = Q + vesicle_center;
	return R;
}

core::Vector
Vesicle::f_transition_f1( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real deriv(f_transition_deriv( conf, resnum, atomnum ) );
	core::Vector r(r_alpha( conf, resnum, atomnum ));

	core::Vector f1_vesicle( f1( xyz, r, deriv));

	if ( !has_pore() ) {
		return f1_vesicle;
	}

	core::Real f_ves( f_vesicle( conf, xyz ) );

	core::Vector f1_p( f1_pore( f_ves, xyz, conf, resnum, atomnum) );
	return f1_vesicle + f1_p;
}

core::Vector
Vesicle::f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real deriv(f_transition_deriv( conf, resnum, atomnum ) );
	core::Vector r(r_alpha( conf, resnum, atomnum ));

	core::Vector f2_vesicle( f2( xyz, r, deriv));

	if ( !has_pore() ) {
		return f2_vesicle;
	}

	core::Real f_ves( f_vesicle( conf, xyz ) );

	core::Vector f2_p( f2_pore( f_ves, xyz, conf, resnum, atomnum) );
	return f2_vesicle + f2_p;
}

//returning string of name of geometry that was created
std::string
Vesicle::geometry_string( ) const {
	std::string geometry_name = "vesicle";
	return geometry_name;
}

//return geometry enum
MP_GEOMETRY_TRANSITION
Vesicle::geometry_enum() const {
	return MP_GEOMETRY_TRANSITION::VESICLE;
}

} // geometry
} // membrane
} // conformation
} // core

