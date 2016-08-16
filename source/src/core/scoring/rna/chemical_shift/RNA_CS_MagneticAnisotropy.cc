// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file   core/scoring/rna/chemical_shift/MagneticAnisotropy.cc
/// @brief
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#include <core/scoring/rna/chemical_shift/RNA_CS_Util.hh>
#include <core/scoring/rna/chemical_shift/RNA_CS_MagneticAnisotropy.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
#include <math.h>
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {

///////////////////////////////////////////////////////////////
///The magnetic_anisotropy contribution of heavy-base atom at src_xyz to the chemqical_shift at CS_data_atom_xyz
Real
delta_magnetic_anisotropy( numeric::xyzVector< core::Real > const & CS_data_atom_xyz,
	numeric::xyzVector< core::Real > const & source_atom_xyz,
	numeric::xyzMatrix< core::Real > const & base_coordinate_matrix,
	RNA_CS_residue_parameters const & source_rsd_CS_params,
	Size const realatomdata_index ) {

	numeric::xyzVector< core::Real > const r_vector = CS_data_atom_xyz - source_atom_xyz;

	Real const r_length = r_vector.length();
	Real const r_length2 = r_length * r_length;
	Real const r_length5 = std::pow( r_length, 5 );

	Real const x_length = dot( r_vector, base_coordinate_matrix.col_x() );
	Real const y_length = dot( r_vector, base_coordinate_matrix.col_y() );
	Real const z_length = dot( r_vector, base_coordinate_matrix.col_z() );

	Real const ma_r = source_rsd_CS_params.magentic_anisotropy_r_coeff(); //diamagnetic component coefficient
	Real const ma_q = source_rsd_CS_params.magentic_anisotropy_q_coeff(); //paramagnetic component coefficient

	Real const coeff_xx = ( ma_r * source_rsd_CS_params.atom_data( realatomdata_index, marx ) ) - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqx ) );
	Real const coeff_yy = ( ma_r * source_rsd_CS_params.atom_data( realatomdata_index, mary ) ) - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqy ) );
	Real const coeff_zz = ( ma_r * source_rsd_CS_params.atom_data( realatomdata_index, marz ) ) - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqz ) );

	Real const coeff_xy =   0.0 - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqw ) );

	//xx component
	Real const termxx = ( ( 3.0 * x_length * x_length ) - r_length2 ) * ( coeff_xx );
	//yy component
	Real const termyy = ( ( 3.0 * y_length * y_length ) - r_length2 ) * ( coeff_yy );
	//zz componenet
	Real const termzz = ( ( 3.0 * z_length * z_length ) - r_length2 ) * ( coeff_zz );

	//xy component
	Real const termxy = ( ( 3.0 * x_length * y_length )             ) * ( coeff_xy );

	Real const termr = -1.0 / ( 3.0 * r_length5 );

	Real const chem_shift_MA = ( termr ) * ( termxx + termyy + termzz + termxy );  //Definition uses opposite sign compare to NUCHEMICS!
	return chem_shift_MA;
}


///////////////////////////////////////////////////////////////
//The gradient of delta_magnetic_anisotropy() with respect to r_vector ( r_vector = CS_data_atom_xyz - source_atom_xyz )
//OK RIGHT NOW CONCERN WITH GETTING THIS FUNCTION RIGHT. OPTIMIZE LATER if NECESSARY!

numeric::xyzVector< core::Real >
get_delta_magnetic_anisotropy_deriv( numeric::xyzVector< core::Real > const & CS_data_atom_xyz,
	numeric::xyzVector< core::Real > const & source_atom_xyz,
	numeric::xyzMatrix< core::Real > const & base_coordinate_matrix,
	RNA_CS_residue_parameters const & source_rsd_CS_params,
	Size const realatomdata_index ) {

	//INCORRECT gradient = (dchem_shift_MA_dx * x_norm) + (dchem_shift_MA_dy * y_norm) + (dchem_shift_MA_dz * z_norm)
	//REASON: base_coordinate_matrix.col_x() and base_coordinate_matrix.col_y() are not entirely orthogonal! (THIS IS DUE TO DEFINITION USED IN NUCNEMICS!)

	//TO ENSURE THAT THE basis are orthogonal:
	//Define yprime_norm=cross(z_norm, x_norm) so that x_norm, yprime_norm and z_norm forms a orthonormal basis of R3. Then:
	//gradient = (dchem_shift_MA_dx * x_norm) + (dchem_shift_MA_dyprime * yprime_norm) + (dchem_shift_MA_dz * z_norm)
	//The crucial thing to remember is to treat y as f(x, y').

	numeric::xyzVector< core::Real > const r_vector = CS_data_atom_xyz - source_atom_xyz;

	Real const r_length = r_vector.length();
	Real const r_length2 = r_length * r_length;
	Real const r_length5 = std::pow( r_length, 5 );

	numeric::xyzVector< core::Real > const x_norm = base_coordinate_matrix.col_x();
	numeric::xyzVector< core::Real > const y_norm = base_coordinate_matrix.col_y();
	numeric::xyzVector< core::Real > const z_norm = base_coordinate_matrix.col_z();
	numeric::xyzVector< core::Real > const yprime_norm = cross( z_norm, x_norm );

	Real const x_length = dot( r_vector, x_norm );
	Real const y_length = dot( r_vector, y_norm );
	Real const z_length = dot( r_vector, z_norm );
	Real const yprime_length = dot( r_vector, yprime_norm );

	//Note that r_length^2= x_length^2 +  yprime_length^2 + z_length^2

	Real const ma_r = source_rsd_CS_params.magentic_anisotropy_r_coeff(); //diamagnetic component coefficient
	Real const ma_q = source_rsd_CS_params.magentic_anisotropy_q_coeff(); //paramagnetic component coefficient

	Real const coeff_xx = ( ma_r * source_rsd_CS_params.atom_data( realatomdata_index, marx ) ) - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqx ) );
	Real const coeff_yy = ( ma_r * source_rsd_CS_params.atom_data( realatomdata_index, mary ) ) - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqy ) );
	Real const coeff_zz = ( ma_r * source_rsd_CS_params.atom_data( realatomdata_index, marz ) ) - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqz ) );

	Real const coeff_xy =   0.0 - ( ma_q * source_rsd_CS_params.atom_data( realatomdata_index, maqw ) );

	//xx component
	Real const termxx = ( ( 3.0 * x_length * x_length ) - r_length2 ) * ( coeff_xx );
	//yy component
	Real const termyy = ( ( 3.0 * y_length * y_length ) - r_length2 ) * ( coeff_yy );
	//zz componenet
	Real const termzz = ( ( 3.0 * z_length * z_length ) - r_length2 ) * ( coeff_zz );

	//xy component
	Real const termxy = ( ( 3.0 * x_length * y_length )             ) * ( coeff_xy );

	Real const termr = -1.0 / ( 3.0 * r_length5 );
	Real const dtermr_dr = ( -5.0 / r_length ) * ( termr );

	//Projection of r_vector onto y-axis gives y= dot(y_norm, r_vector) =  dot(y_norm, yprime_norm) * y_prime + dot(y_norm, xprime_norm) * x
	Real const dy_dx  = dot( y_norm, x_norm );
	Real const dy_dyprime = dot( y_norm, yprime_norm );

	/////////////////////////////////
	//dchem_shift_MA_dx
	/////////////////////////////////

	//dchem_shift_MA_dx = ( (dtermr_dx) * (termxx + termyy + termzz + termxy ) ) + ( (termr) * (dtermxx_dx + dtermyy_dx + dtermzz_dx + dtermxy_dx) )   //Chain rule.

	Real const dr_dx = ( x_length / r_length );
	Real const dtermr_dx  = dr_dx * dtermr_dr;
	Real const dtermxx_dx = ( ( 6.0 * x_length                 )  - ( 2.0 * x_length ) ) * ( coeff_xx );
	Real const dtermyy_dx = ( ( 6.0 * y_length * dy_dx             ) - ( 2.0 * x_length ) ) * ( coeff_yy );
	Real const dtermzz_dx = ( ( 0.0                            ) - ( 2.0 * x_length ) ) * ( coeff_zz );
	Real const dtermxy_dx = ( ( ( 3.0 * y_length ) + ( 3.0 * x_length * dy_dx ) ) - ( 0.0            ) ) * ( coeff_xy );

	Real const dchem_shift_MA_dx = ( ( dtermr_dx ) * ( termxx + termyy + termzz + termxy ) ) + ( ( termr ) * ( dtermxx_dx + dtermyy_dx + dtermzz_dx + dtermxy_dx ) );

	/////////////////////////////////
	//dchem_shift_MA_dyprime
	/////////////////////////////////

	//dchem_shift_MA_dyprime = ( (dtermr_dyprime) * (termxx + termyy + termzz + termxy ) )
	//             + ( (termr) * (dtermxx_dyprime + dtermyy_dyprime + dtermzz_dyprime + dtermxy_dyprime) )   //Chain rule.

	Real const dr_dyprime = ( yprime_length / r_length );
	Real const dtermr_dyprime  = dr_dyprime * dtermr_dr;
	Real const dtermxx_dyprime = ( ( 0.0                    ) - ( 2.0 * yprime_length ) ) * ( coeff_xx );
	Real const dtermyy_dyprime = ( ( 6.0 * y_length * dy_dyprime ) - ( 2.0 * yprime_length ) ) * ( coeff_yy );
	Real const dtermzz_dyprime = ( ( 0.0                   ) - ( 2.0 * yprime_length ) ) * ( coeff_zz );
	Real const dtermxy_dyprime = ( ( 3.0 * x_length * dy_dyprime ) - ( 0.0              ) ) * ( coeff_xy );

	Real const dchem_shift_MA_dyprime = ( ( dtermr_dyprime ) * ( termxx + termyy + termzz + termxy ) ) + ( ( termr ) * ( dtermxx_dyprime + dtermyy_dyprime + dtermzz_dyprime + dtermxy_dyprime ) );

	/////////////////////////////////
	//dchem_shift_MA_dz
	/////////////////////////////////

	//dchem_shift_MA_dz = ( (dtermr_dz) * (termxx + termyy + termzz + termxy ) ) + ( (termr) * (dtermxx_dz + dtermyy_dz + dtermzz_dz + dtermxy_dz) )   //Chain rule.

	Real const dr_dz = ( z_length / r_length );
	Real const dtermr_dz  = dr_dz * dtermr_dr;
	Real const dtermxx_dz = ( ( 0.0            ) - ( 2.0 * z_length ) ) * ( coeff_xx );
	Real const dtermyy_dz = ( ( 0.0            ) - ( 2.0 * z_length ) ) * ( coeff_yy );
	Real const dtermzz_dz = ( ( 6.0 * z_length ) - ( 2.0 * z_length ) ) * ( coeff_zz );
	Real const dtermxy_dz = ( ( 0.0            ) - ( 0.0            ) ) * ( coeff_xy );

	Real const dchem_shift_MA_dz = ( ( dtermr_dz ) * ( termxx + termyy + termzz + termxy ) ) + ( ( termr ) * ( dtermxx_dz + dtermyy_dz + dtermzz_dz + dtermxy_dz ) );

	/////////////////////////////////////////////////////////////////////

	numeric::xyzVector< core::Real > const analytical_gradient = ( dchem_shift_MA_dx * x_norm ) + ( dchem_shift_MA_dyprime * yprime_norm ) + ( dchem_shift_MA_dz * z_norm );

	/////////////////////////////////////////////////////////////////////

	return analytical_gradient;
}

///////////////////////////////////////////////////////////////
///The magnetic_anisotropy contribution of source_rsd to the chemical_shift at atom_xyz

Real
magnetic_anisotropy_effect( numeric::xyzVector< core::Real > const & atom_xyz, conformation::Residue const & source_rsd, RNA_CS_residue_parameters const & source_rsd_CS_params ){

	if ( source_rsd.aa() != source_rsd_CS_params.aa() ) utility_exit_with_message( "rsd.aa() != source_rsd_CS_params.aa()!" );

	Real ma_effect = 0.0;
	Size const maxatoms = source_rsd_CS_params.get_atomnames_size();

	numeric::xyzMatrix< core::Real > const base_coordinate_matrix = get_rna_base_coordinate_system_from_CS_params( source_rsd, source_rsd_CS_params );

	for ( Size realatomdata_index = 1; realatomdata_index < maxatoms; realatomdata_index++ ) {
		if ( dround( source_rsd_CS_params.atom_data( realatomdata_index, maca ) ) != 1 ) continue;
		Size const atom_index = source_rsd.atom_index( source_rsd_CS_params.get_atomname( realatomdata_index ) );
		ma_effect += delta_magnetic_anisotropy( atom_xyz, source_rsd.xyz( atom_index ), base_coordinate_matrix, source_rsd_CS_params, realatomdata_index );
	}

	return ma_effect;
}

///////////////////////////////////////////////////////////////


} //chemical_shift
} //rna
} //scoring
} //core


