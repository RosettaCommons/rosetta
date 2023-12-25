// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/membrane_geometry/Slab.cc
/// @brief Data describing the parameters of a flat membrane
///
/// @details Slab class contains the parameters of a flat membrane and
///  the function to calculate the transition from the hydrophobic
///  environment of the membrane to a hydrophilic environment.
///
/// @note This object is a member of Conformation and should only be accessed using
///            pose.conformation().membrane_geometry().
///
/// @author Hope Woods (hope.woods@vanderbilt.edu)
/// @note Adding the electrostatic part for slab geometry on Nov 30th, 2023.
/// @author Rituparna Samanta (rituparna@utexas.edu)
// Unit Headers
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/membrane_geometry/Slab.hh>

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

// C++ Headers
#include <string>

static basic::Tracer TR( "core.conformation.membrane.membrane_geometry.Slab" );

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

Slab::Slab(
	core::Real steepness
) :
	MembraneGeometry( steepness )
{}

Slab::Slab(
	core::Real steepness,
	core::Real thickness
) :
	MembraneGeometry( steepness, thickness )
{}


Slab::Slab(
	core::Real steepness,
	core::Real thickness,
	AqueousPoreParametersOP aqueous_pore
) :
	MembraneGeometry( steepness, thickness, aqueous_pore )
{}

/// @brief Destructor
Slab::~Slab() {}

MembraneGeometryOP Slab::clone() const {
	return SlabOP( new Slab( *this ) );
}

/// @brief Generate a string representation of information represented by Slab
void
Slab::show() const {
	show( std::cout );
}

void
Slab::show( std::ostream & output ) const {
	output << "MembraneGeometry: Information about the geometry of the membrane or membrane mimetic" << std::endl;
	//output << "Membrane Thickness: " << thickness_ << std::endl;
	//output << "Membrane Steepness: " << steepness_ << std::endl;
}

//Slab transition function and helper functions
//xyz is the coordinates in space of the atom of interest
//n is steepness of hydrophobic -> hydrophillic transition (default = 15)
//returns the value of the transition function for membrane score functions
core::Real
Slab::f_transition( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );

	core::Real f_thk( f_thickness( conf, xyz.z()) );
	return f_hydration( f_thk, xyz );

}

core::Real
Slab::f_transition_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {

	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );

	core::Real f_thk_deriv_dz( f_thickness_deriv( conf, xyz.z() ) );

	if ( !has_pore() ) {
		return f_thk_deriv_dz;
	} else {
		return f_hydration_deriv_dz( xyz, f_thk_deriv_dz );
	}

}

//returns electrostatic field due to the slab geometry
//xyz is the coordinates in space of the atom of interest
core::Real
Slab::fa_elec_lipid( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	//TR << "lipid name: " << conf.membrane_info()->implicit_lipids()->lipid_composition_name() << std::endl ;
	core::Real f_elec_lipid( conf.membrane_info()->implicit_lipids()->f_elec_field( xyz.z() ) );
	return f_elec_lipid;

}

core::Real
Slab::fa_elec_lipid_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {

	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );

	//TR << "lipid name: " << conf.membrane_info()->implicit_lipids()->lipid_composition_name() << std::endl ;
	core::Real f_elec_lipid_deriv( conf.membrane_info()->implicit_lipids()->f_elec_field_gradient( xyz.z() ) );
	return f_elec_lipid_deriv;

}

core::Vector
Slab::r_alpha( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Real z_depth = conf.membrane_info()->atom_z_position( conf, resnum, atomnum );
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	const core::Vector normal = conf.membrane_info()->membrane_normal( conf );
	core::Vector const & xyz( conf.residue( resnum ).atom( atomnum ).xyz() );
	core::Vector proj_i = mem_cen + z_depth * normal;
	core::Vector i_ip = proj_i - xyz;
	return ( mem_cen - i_ip );
}

// Extracellular
// .............................................................................
//                                             * <- atom(m)
//                                             |
//                                             |
// -----Membrane Center--------Origin->*-------* <-f_alpha------------------
//
//
//
// .............................................................................
// Intracelluar

core::Vector
Slab::f_transition_f1( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real deriv(f_transition_deriv( conf, resnum, atomnum ) );
	core::Vector r(r_alpha( conf, resnum, atomnum ));
	core::Vector f1_z( f1( xyz, r, deriv));

	if ( !has_pore() ) {
		return f1_z;
	}

	core::Real f_thk( f_thickness( conf, xyz.z()) );
	core::Vector f1_p( f1_pore( f_thk, xyz, conf, resnum, atomnum) );

	return f1_z + f1_p;
}

core::Vector
Slab::f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real deriv(f_transition_deriv( conf, resnum, atomnum ) );
	core::Vector r(r_alpha( conf, resnum, atomnum ));
	core::Vector f2_z( f2( xyz, r, deriv));
	if ( !has_pore() ) {
		return f2_z;
	}
	core::Real f_thk( f_thickness( conf, xyz.z()) );
	core::Vector f2_p( f2_pore( f_thk, xyz, conf, resnum, atomnum) );

	return f2_z + f2_p;
}


//returning string of name of geometry that was created
std::string
Slab::geometry_string( ) const {
	std::string geometry_name = "slab";
	return geometry_name;
}

//return geometry enum
MP_GEOMETRY_TRANSITION
Slab::geometry_enum() const {
	return MP_GEOMETRY_TRANSITION::SLAB;
}

} // geometry
} // membrane
} // conformation
} // core

