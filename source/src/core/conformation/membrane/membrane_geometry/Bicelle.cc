// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/membrane_geometry/Bicelle.cc
/// @brief Data describing the parameters of the bicelle
///
/// @details Bicelle class contains the parameters of the bicelle and
///  the function to calculate the transition from the hydrophobic
///  environment of the bicelle to a hydrophilic environment.
///
/// @note This object is a member of Conformation and should only be accessed using
///            pose.conformation().membrane_geometry().
/// @author Hope Woods (hope.woods@vanderbilt.edu)

// Unit Headers
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/membrane_geometry/Bicelle.hh>

// Project Headers
#include <core/conformation/membrane/MembraneParams.hh>

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

static basic::Tracer TR( "core.conformation.membrane.membrane_geometry.Bicelle" );

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


//To initialize an instance of Bicelle, either the conformation and membrane pos
//must be provided or the Bicelle inner radius.
//The conformation and membrane position are needed to estimate the diamter of the protein
//at the center of the membrane.
Bicelle::Bicelle(
	core::Real steepness,
	Conformation const & conf,
	core::Size membrane_pos
) :
	MembraneGeometry( steepness )
{
	core::Real protein_diameter = protein_slice_diameter_at_mem_cen( conf, membrane_pos );
	set_protein_slice_diameter( protein_diameter);
}

Bicelle::Bicelle(
	core::Real steepness,
	core::Real thickness,
	Conformation const & conf,
	core::Size membrane_pos

) :
	MembraneGeometry( steepness, thickness )
{

	core::Real protein_diameter = protein_slice_diameter_at_mem_cen( conf, membrane_pos );
	set_protein_slice_diameter( protein_diameter);
}

Bicelle::Bicelle(
	core::Real steepness,
	core::Real thickness,
	core::Real bicelle_inner_radius
) :
	MembraneGeometry( steepness, thickness )
{
	set_inner_radius( bicelle_inner_radius );
	update_outer_radius();
}

/// @brief Destructor
Bicelle::~Bicelle() {}


BicelleOP
Bicelle::clone() const {
	return BicelleOP( new Bicelle( *this ) );
}

/// @brief Generate a string representation of information represented by Bicelle
void
Bicelle::show() const {
	show( std::cout );
}

void
Bicelle::show( std::ostream & output ) const {
	output << "MembraneGeometry: Information about the geometry of the membrane or membrane mimetic" << std::endl;
	output << "Bicelle Inner Radius: " << bicelle_inner_radius_ << std::endl;
}

core::Real
Bicelle::bicelle_inner_radius() const{
	return bicelle_inner_radius_;
}


core::Real
Bicelle::protein_slice_diameter_at_mem_cen( Conformation const & conf, core::Size membrane_pos ) const {

	core::Size num_res = conf.size();
	//array to store residue positions that are within 3A of membrane center
	utility::vector1<core::Vector> res_xyz;
	//array to store distances between CA atoms
	utility::vector1<core::Real> dist_array;
	//grab membrane_normal vector
	core::Vector mem_normal = conf.residue( membrane_pos ).atom( membrane::normal ).xyz();
	mem_normal.normalize();
	//grab membrane_thickness vector
	core::Vector mem_thick = conf.residue( membrane_pos ).atom( membrane::thickness ).xyz();
	mem_thick.normalize();
	//binormal_vector
	core::Vector mem_binormal = mem_normal.cross( mem_thick );
	//grab membrane center
	core::Vector mem_center = conf.residue( membrane_pos ).atom( membrane::center ).xyz();
	//loop through residues
	for ( core::Size i = 1; i <= num_res; i++ ) {
		//subtract membrane center from xyz coord
		core::Vector diff = conf.residue( i ).atom( 2 ).xyz() - mem_center;
		//dot product of difference and membrane normal
		core::Real dist_from_center_i = dot( diff, mem_normal);

		//if CA within 3 A of membrane center continue
		if ( std::abs(dist_from_center_i) <= 3 ) {
			res_xyz.push_back( conf.residue(i).atom(2).xyz() );
		}

		//loop through residues within 3 A of membrane center
		for ( core::Size ai = 1; ai <= res_xyz.size(); ai++ ) {
			//loop through list again
			for ( core::Size aj = 1; aj <= res_xyz.size(); aj++ ) {
				//skip if same resiue
				if ( aj != ai ) {
					core::Vector i_xyz = res_xyz[ai];
					core::Vector j_xyz = res_xyz[aj];
					//calcuate distance between CA atoms in xy dimensions
					core::Real dist = pow( pow(corrected_coordinate( i_xyz, mem_thick ) - corrected_coordinate( j_xyz, mem_thick ), 2) + pow(corrected_coordinate( i_xyz, mem_binormal ) - corrected_coordinate( j_xyz, mem_binormal ), 2), 0.5);
					//store in array
					dist_array.push_back( dist );
				}
			}
		}
	}
	//find 95 percentile of distances
	std::sort(dist_array.begin(), dist_array.end());
	//multiply k percent by the total number of values,n
	core::Size index = 0.95 * dist_array.size();
	core::Real protein_core = dist_array[index];
	std::cout << "protein_core: " << protein_core << std::endl;
	return protein_core;
}


//Sets protein_slice_diameter_, should be calculated by protein_slice_diameter_at_mem_cen
//calls update_radii() to update bicelle_inner_radius and bicelle_outer_radius after protein_slice_diameter_ is set.
void
Bicelle::set_protein_slice_diameter( core::Real diameter ) {
	protein_slice_diameter_ = diameter;
	TR << "setting_protein_slice_diameter for Bicelle: " << protein_slice_diameter_ << std::endl;
	update_radii(); //radius depends on protein_slice_diameter so if set_protein_slice_diameter is called update_radii should also run
}

void
Bicelle::set_inner_radius( core::Real inner_r ) {
	if ( inner_r < 0 ) {
		TR.Fatal << "Cannot set Bicelle inner radius as a negative number." << std::endl;
	} else if ( protein_slice_diameter_ == 0.0 ) {
		TR.Warning << "Protein slice diameter was not set or unable to find residues at membrane center." << std::endl;
		TR.Warning << "Cannot be sure bicelle is a reasonable size for the protein." << std::endl;
	}
	bicelle_inner_radius_ = inner_r;
}

void
Bicelle::set_outer_radius( core::Real outer_r ) {
	if  ( outer_r < bicelle_inner_radius_ ) {
		TR.Fatal << "Bicelle outer radius must be greater than the set inner radius." << std::endl;
	}
	bicelle_outer_radius_ = outer_r;
	update_edge_steepness(); //if outer_radius is changed, edge_steepness needs to be updated
}

void
Bicelle::set_bicelle_edge_steepness( core::Real edge_steepness ) {
	bicelle_edge_steepness_ = edge_steepness;
}

core::Real
Bicelle::protein_slice_diameter() const {
	return protein_slice_diameter_;
}



void Bicelle::update_inner_radius() {

	core::Real recommended_micelle_inner_radius = 3*(protein_slice_diameter_/2);
	//did the user set the radius?
	if ( basic::options::option[ basic::options::OptionKeys::mp::geo::bicelle_radius ].user() ) {
		set_inner_radius( basic::options::option[ basic::options::OptionKeys::mp::geo::bicelle_radius ]());
		TR << "setting inner radius from command line option:  " << bicelle_inner_radius_ << std::endl;
		if ( protein_slice_diameter_ >= (bicelle_inner_radius_*2) ) {
			TR.Warning << "Protein core at the membrane center is larger than the radius of your bicelle/micelle." << std::endl;
		}
		if ( bicelle_inner_radius_ < recommended_micelle_inner_radius ) {
			TR.Warning << "Consider setting a larger radius. The set radius is smaller than the recommended radius of " << recommended_micelle_inner_radius << std::endl;
		}
	} else {
		TR << "protein_slice_diameter: " << protein_slice_diameter_ << std::endl;
		set_inner_radius( recommended_micelle_inner_radius );
		TR << "setting inner radius as: " << bicelle_inner_radius_ << std::endl;
	}
}

//set outer radius based on inner radius and thickness
void Bicelle::update_outer_radius() {
	core::Real mem_thickness = membrane_thickness();
	core::Real inner_radius = bicelle_inner_radius_;
	set_outer_radius( inner_radius + mem_thickness );
}

core::Real
Bicelle::bicelle_outer_radius() const{
	return bicelle_outer_radius_;
}

void Bicelle::update_radii() {
	update_inner_radius();
	update_outer_radius();
}


// For the bicelle there are two steepnesses. One in the z-direction (membrane normal)
// and one in the x-y direction. The membrane_steepness is the one in the z-direction;
// the steepness in xy direction should be set according to the dimension in the xy
// plane which includes the bicelle radius.
//set bicelle_edge_steepness
void Bicelle::update_edge_steepness() {
	core::Real mem_thickness = membrane_thickness();
	core::Real mem_steepness = membrane_steepness();
	core::Real outer_radius = bicelle_outer_radius_;
	core::Real a = (outer_radius + mem_thickness*2)/(mem_thickness*2);
	set_bicelle_edge_steepness( mem_steepness*a );
}

//get bicelle_edge_steepness
core::Real Bicelle::bicelle_edge_steepness() const {
	return bicelle_edge_steepness_;
}


//Bicelle transition function and helper functions

//h_bicelle is transition functino for bicelle edge
//xyz is the coordinates in space of the atom of interest
//n is steepness of hydrophobic -> hydrophillic transition
core::Real
Bicelle::h_bicelle( core::Vector xyz, const core::Vector mem_cen ) const {
	core::Real r = xyz.distance(mem_cen); //mem_cen in this case is a point at the center of the bicelle
	core::Real n = bicelle_edge_steepness_;
	core::Real rp = r/bicelle_outer_radius_;
	core::Real h_r = pow(rp, n)/(1+ pow(rp, n));
	return h_r;
}


core::Real
Bicelle::h_bicelle_deriv_wrt_r( core::Vector xyz, const core::Vector mem_cen ) const {
	core::Real r = xyz.distance(mem_cen);
	core::Real rp = r/bicelle_outer_radius_;
	core::Real n = bicelle_edge_steepness_;
	core::Real numerator = n*pow(rp, n-1);
	core::Real denominator = bicelle_outer_radius_*pow( 1+pow(rp, n), 2);
	return (numerator/denominator);
}

//combine h_bicelle and f_imm1
//f is the value of the transition function, while f_z and h_r are components
//that make up f.
core::Real
Bicelle::f_bicelle( core::Vector xyz, core::Real z_depth, const core::Vector mem_cen ) const {
	core::Real h_r = h_bicelle( xyz, mem_cen );
	core::Real f_z = f_imm1( z_depth );
	core::Real f = f_z + h_r - ( f_z*h_r );
	return f;
}

//Same as above where f_z and h_r are components of the transition function:
//f = f_z + h_r - f_z*h_r
core::Real
Bicelle::f_bicelle_deriv( core::Vector xyz, core::Real z_depth, const core::Vector mem_cen ) const {
	core::Real h_r = h_bicelle( xyz, mem_cen );
	core::Real f_z = f_imm1( z_depth );
	core::Real h_deriv = h_bicelle_deriv_wrt_r( xyz, mem_cen );
	core::Real f_z_deriv = f_imm1_deriv( z_depth );

	core::Real bicelle_deriv = f_z_deriv + h_deriv - (f_z_deriv*h_r + f_z*h_deriv);

	return bicelle_deriv;
}

core::Real
Bicelle::f_transition( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real z_depth = conf.membrane_info()->atom_z_position( conf, resnum, atomnum );
	return f_bicelle( xyz, z_depth, mem_cen );
}

core::Real
Bicelle::f_transition_deriv_m( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Real z_depth = conf.membrane_info()->atom_z_position( conf, resnum, atomnum );
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real df_dz = f_imm1_deriv( z_depth );
	core::Real h = h_bicelle( xyz, mem_cen );
	return df_dz*( 1 - h );
}

core::Real
Bicelle::f_transition_deriv( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	core::Real z_depth = conf.membrane_info()->atom_z_position( conf, resnum, atomnum );
	core::Real dh_dr = h_bicelle_deriv_wrt_r( xyz, mem_cen );
	core::Real f = f_imm1( z_depth );
	return dh_dr*( 1 - f );
}

core::Vector
Bicelle::r_alpha_m( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Real z_depth = conf.membrane_info()->atom_z_position( conf, resnum, atomnum );
	const core::Vector mem_cen = conf.membrane_info()->membrane_center( conf );
	const core::Vector normal = conf.membrane_info()->membrane_normal( conf );
	core::Vector const & xyz( conf.residue( resnum ).atom( atomnum ).xyz() );
	core::Vector proj_i = mem_cen + z_depth * normal;
	core::Vector i_ip = proj_i - xyz;
	return ( mem_cen - i_ip );
}

core::Vector
Bicelle::r_alpha( Conformation const & conf, core::Size, core::Size ) const {
	return conf.membrane_info()->membrane_center( conf );
}


core::Vector
Bicelle::f_transition_f1( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	//f1 for imm1, z dependence
	core::Real deriv_m( f_transition_deriv_m( conf, resnum, atomnum ));
	core::Vector r_m(r_alpha_m( conf, resnum, atomnum ));
	core::Vector f1_m = f1( xyz, r_m, deriv_m);

	//f1 for bicelle, xyz dependence
	core::Real deriv_bic( f_transition_deriv( conf, resnum, atomnum ));
	core::Vector r_bic(r_alpha( conf, resnum, atomnum ));
	core::Vector f1_bicelle = f1( xyz, r_bic, deriv_bic);
	return f1_m + f1_bicelle;
}

core::Vector
Bicelle::f_transition_f2( Conformation const & conf, core::Size resnum, core::Size atomnum ) const {
	core::Vector const & xyz( corrected_xyz( conf, resnum, atomnum) );
	//f2 for imm1, z dependence
	core::Real deriv_m( f_transition_deriv_m( conf, resnum, atomnum ));
	core::Vector r_m(r_alpha_m( conf, resnum, atomnum ));
	core::Vector f2_m = f2( xyz, r_m, deriv_m);

	//f2 for bicelle, xyz dependence
	core::Real deriv_bic( f_transition_deriv( conf, resnum, atomnum ));
	core::Vector r_bic(r_alpha( conf, resnum, atomnum ));
	core::Vector f2_bicelle = f2( xyz, r_bic, deriv_bic);
	return f2_m + f2_bicelle;
}


//returning string of name of geometry that was created
std::string
Bicelle::geometry_string( ) const {
	std::string geometry_name = "bicelle";
	return geometry_name;
}

//return geometry enum
MP_GEOMETRY_TRANSITION
Bicelle::geometry_enum() const {
	return MP_GEOMETRY_TRANSITION::BICELLE;
}


} // geometry
} // membrane
} // conformation
} // core

