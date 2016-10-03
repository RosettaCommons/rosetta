// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file   core/scoring/rna/chemical_shift/RNA_CS_RingCurrent.cc
/// @brief
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#include <core/scoring/rna/chemical_shift/RNA_CS_Util.hh>
#include <core/scoring/rna/chemical_shift/RNA_CS_RingCurrent.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {

///////////////////////////////////////////////////////////////
numeric::xyzVector< core::Real >
ring_pos( conformation::Residue const & rsd, RNA_CS_residue_parameters const & rna_cs_rsd_params, Size const ring_ID ) //NOTE: ONLY INCLUDE HEAVY ATOMS (no hydrogens!)
{
	numeric::xyzVector< core::Real > ring_center( 0.0, 0.0, 0.0 );

	Size atom_count = 0;

	if ( rsd.aa() != rna_cs_rsd_params.aa() ) utility_exit_with_message( "rsd.aa() != rna_cs_rsd_params.aa()!" );

	Size const maxatoms = rna_cs_rsd_params.get_atomnames_size();

	for ( Size count = 1; count < maxatoms; count++ ) {
		if ( ring_ID == 1 ) {
			if ( dround( rna_cs_rsd_params.atom_data( count, rcl1 ) != 1 ) ) continue;
		} else if ( ring_ID == 2 ) {
			if ( dround( rna_cs_rsd_params.atom_data( count, rcl2 ) != 1 ) ) continue;
		} else {
			utility_exit_with_message( "( ring_ID != 1 ) && ( ring_ID != 2 ), ring_ID = ( " + ObjexxFCL::string_of( ring_ID ) + " )!" );
		}

		Size const atom_index = rsd.atom_index( rna_cs_rsd_params.get_atomname( count ) );

		ring_center += rsd.xyz( atom_index );
		atom_count++;
	}

	if ( atom_count == 0 ) utility_exit_with_message( "atom_count == 0!" );
	ring_center = ring_center/atom_count;
	return ring_center;
}

///////////////////////////////////////////////////////////////
/* derivative of ellintk. expects 0 <= m < 1 as input.
*/
static
Real dellintk_dm( Real const m )
{
	//Real const ak0 = 1.38629436112;
	Real const ak1 = 0.09666344259;
	Real const ak2 = 0.03590092383;
	Real const ak3 = 0.03742563713;
	Real const ak4 = 0.01451196212;
	Real const bk0 = 0.5          ;
	Real const bk1 = 0.12498593597;
	Real const bk2 = 0.06880248576;
	Real const bk3 = 0.03328355346;
	Real const bk4 = 0.00441787012;
	Real const mr = 1.0 - m;
	Real const mr2 = mr * mr;
	Real const mr3 = mr * mr * mr;
	Real const mr4 = mr * mr * mr * mr;

	Real const bd = bk0 + ( bk1 * mr ) + ( bk2 * mr2 ) + ( bk3 * mr3 ) + ( bk4 * mr4 );
	Real const cd = std::log ( 1.0 / mr );

	//DEFINITION: ellintk        =     ad + (bd * cd)
	//            dellinte_dm    =     dad_dm + (bd * dcd_dm + dbd_dm * cd)

	Real const dad_dm = -1.0 * ( ( ak1 + ( 2 * ak2 * mr ) + ( 3 * ak3 * mr2 ) + ( 4 * ak4 * mr3 ) ) );
	Real const dbd_dm = -1.0 * ( ( bk1 + ( 2 * bk2 * mr ) + ( 3 * bk3 * mr2 ) + ( 4 * bk4 * mr3 ) ) );
	Real const dcd_dm = + 1.0 * ( 1.0 / mr  );

	return ( dad_dm + ( bd * dcd_dm + dbd_dm * cd ) );
}

///////////////////////////////////////////////////////////////
/* derivative of ellinte. expects 0 <= m < 1 as input.
*/
static
Real dellinte_dm( Real const m )
{
	//Real const ae0 = 1.0          ;
	Real const ae1 = 0.44325141463;
	Real const ae2 = 0.06260601220;
	Real const ae3 = 0.04757383546;
	Real const ae4 = 0.01736506451;
	Real const be0 = 0.0          ;
	Real const be1 = 0.24998368310;
	Real const be2 = 0.09200180037;
	Real const be3 = 0.04069697526;
	Real const be4 = 0.00526449639;
	Real const mr = 1.0 - m;
	Real const mr2 = mr * mr;
	Real const mr3 = mr * mr * mr;
	Real const mr4 = mr * mr * mr * mr;

	Real const bd = be0 + ( be1 * mr ) + ( be2 * mr2 ) + ( be3 * mr3 ) + ( be4 * mr4 );
	Real const cd = std::log( 1.0 / mr );

	//DEFINITION: ellinte        =     ad + (bd * cd)
	//            dellinte_dm    =     dad_dm + (bd * dcd_dm + dbd_dm * cd)

	Real const dad_dm = -1.0 * ( ( ae1 + ( 2 * ae2 * mr ) + ( 3 * ae3 * mr2 ) + ( 4 * ae4 * mr3 ) ) );
	Real const dbd_dm = -1.0 * ( ( be1 + ( 2 * be2 * mr ) + ( 3 * be3 * mr2 ) + ( 4 * be4 * mr3 ) ) );
	Real const dcd_dm = + 1.0 * ( 1.0 / mr  );

	return ( dad_dm + ( bd * dcd_dm + dbd_dm * cd ) );
}

///////////////////////////////////////////////////////////////
/* ellintk, calculates the K factor of the elliptical integral
expects 0 <= m < 1 as input. This function was difined from formula 17.3.34
Abramowitz and Segun, Handbook of Mathematical Functions,
Dover publications 1965
*/
static
Real ellintk( Real const m ) //ellintk( m = 0 ) = PI/2
{
	Real const ak0 = 1.38629436112;
	Real const ak1 = 0.09666344259;
	Real const ak2 = 0.03590092383;
	Real const ak3 = 0.03742563713;
	Real const ak4 = 0.01451196212;
	Real const bk0 = 0.5          ;
	Real const bk1 = 0.12498593597;
	Real const bk2 = 0.06880248576;
	Real const bk3 = 0.03328355346;
	Real const bk4 = 0.00441787012;
	Real const mr = 1.0 - m;
	Real const mr2 = mr * mr;
	Real const mr3 = mr * mr * mr;
	Real const mr4 = mr * mr * mr * mr;

	Real const ad = ak0 + ( ak1 * mr ) + ( ak2 * mr2 ) + ( ak3 * mr3 ) + ( ak4 * mr4 );
	Real const bd = bk0 + ( bk1 * mr ) + ( bk2 * mr2 ) + ( bk3 * mr3 ) + ( bk4 * mr4 );
	Real const cd = std::log ( 1.0 / mr );

	return ( ad + ( bd * cd ) );
}


///////////////////////////////////////////////////////////////
/* ellinte, calculates the E factor of the elliptical integral
expects 0 <= m < 1 as input. This function was defined from formula 17.3.36
Abramowitz and Segun, Handbook of Mathematical Functions,
Dover publications 1965
*/
static
Real ellinte( Real const m ) //ellinte( m = 0 ) = PI/2
{
	Real const ae0 = 1.0          ;
	Real const ae1 = 0.44325141463;
	Real const ae2 = 0.06260601220;
	Real const ae3 = 0.04757383546;
	Real const ae4 = 0.01736506451;
	Real const be0 = 0.0          ;
	Real const be1 = 0.24998368310;
	Real const be2 = 0.09200180037;
	Real const be3 = 0.04069697526;
	Real const be4 = 0.00526449639;
	Real const mr = 1.0 - m;
	Real const mr2 = mr * mr;
	Real const mr3 = mr * mr * mr;
	Real const mr4 = mr * mr * mr * mr;

	Real const ad = ae0 + ( ae1 * mr ) + ( ae2 * mr2 ) + ( ae3 * mr3 ) + ( ae4 * mr4 );
	Real const bd = be0 + ( be1 * mr ) + ( be2 * mr2 ) + ( be3 * mr3 ) + ( be4 * mr4 );
	Real const cd = std::log ( 1.0 / mr );

	return ( ad + ( bd * cd ) );
}

///////////////////////////////////////////////////////////////
/*calculates the cylindrical coordinates for atom_xyz with the origin located ring_center and z - direction equal to base_z_axis  */
void
get_rho_and_z( numeric::xyzVector< core::Real > const & atom_xyz, numeric::xyzVector< core::Real > const &  ring_center, numeric::xyzVector< core::Real > const & base_z_axis, Real & rho, Real & z )
{
	numeric::xyzVector< core::Real > const r_vector = atom_xyz - ring_center;

	Real const r_length = r_vector.length();

	if ( std::abs( base_z_axis.length()  - 1.0 )   > 0.00001 ) utility_exit_with_message( "std::abs( base_z_axis.length()  - 1.0 )   > 0.00001 !!" );

	Real const dot_product = dot( r_vector, base_z_axis ) / ( r_length ); //i.e cos( angle )
	Real const angle = acos( dot_product );

	rho = sin ( angle ) * r_length;
	z   = dot_product * r_length;
}

///////////////////////////////////////////////////////////////
Real
delta_ring_current_term( Real const rho, Real const z, Real const ring_radius, Real const ring_z_offset )
{
	Real const d = rho / ring_radius;
	Real const hp = ( z - ring_z_offset ) / ring_radius;
	Real const hp2 = hp * hp;

	Real const rp = ( ( 1.0 + d ) * ( 1.0 + d ) ) + hp2;
	Real const sp = ( ( 1.0 - d ) * ( 1.0 - d ) ) + hp2;
	Real const tp = ( 1.0 - ( d * d ) - hp2 );

	Real const mp = ( 4.0 * d ) / rp;

	Real const Fp = ( 2.0 / std::sqrt ( rp ) );
	Real const Gp = ( tp / sp );
	Real const Hp = ( ellintk( mp ) + ( Gp * ellinte( mp ) ) );

	Real const term = ( Fp ) * ( Hp );

	return term;
}

///////////////////////////////////////////////////////////////
void
get_ring_current_term_derivatives( Real const rho, Real const z, Real const ring_radius, Real const ring_z_offset, Real & dterm_drho, Real & dterm_dz )
{
	Real const d = rho / ring_radius;
	Real const hp = ( z - ring_z_offset ) / ring_radius;
	Real const hp2 = hp * hp;

	Real const rp = ( ( 1.0 + d ) * ( 1.0 + d ) ) + hp2;
	Real const sp = ( ( 1.0 - d ) * ( 1.0 - d ) ) + hp2;
	Real const tp = ( 1.0 - ( d * d ) - hp2 );

	Real const rp2 = rp * rp;
	Real const rp3 = rp * rp * rp;

	Real const sp2 = sp * sp;

	Real const mp = ( 4.0 * d ) / rp;

	Real const Fp = ( 2.0 / std::sqrt ( rp ) );
	Real const Gp = ( tp / sp );
	Real const Hp = ( ellintk( mp ) + ( Gp * ellinte( mp ) ) );

	//Definition: term = ( Fp ) * ( Hp );

	/////////////////////////////////
	//dterm_drho
	/////////////////////////////////

	//dpterm_drho = (dFp_drho * Hp) + (Fp * dHp_drho) //Chain rule.

	Real const dFp_drho = ( 1.0 / ring_radius ) * (  - 2.0 * ( 1.0 + d ) ) * ( 1 / ( sqrt ( rp3 ) ) );

	//EQ: dHp_drho = d_ellintk_drho + (dGp_drho * ellinte(mp) ) + ( Gp * d_ellinte_drho )
	//EQ: dHp_drho = (d_ellintk_dmp * dmp_drho )  + (dGp_drho * ellinte(mp) ) + ( Gp * dellinte_dmp * dmp_drho )

	Real const dmp_drho = ( 1.0 / ring_radius ) * ( ( rp * 4.0 ) - ( ( 4.0 * d ) * ( 2.0 * ( 1.0 + d ) ) ) ) * ( 1.0 / ( rp2 ) ); //Quotient Rule.

	Real const dGp_drho = ( 1.0 / ring_radius ) * ( ( sp * ( -2.0 * d ) ) - ( tp * (  - 2.0 * ( 1.0 - d ) ) ) ) * ( 1.0 / ( sp2 ) ); //Quotient Rule.

	Real const dHp_drho = ( dellintk_dm( mp ) * dmp_drho )  + ( dGp_drho * ellinte( mp ) ) + ( Gp * dellinte_dm( mp ) * dmp_drho );

	dterm_drho = ( dFp_drho * Hp ) + ( Fp * dHp_drho );

	/////////////////////////////////
	//dterm_dz
	/////////////////////////////////

	//dpterm_dz = (dFp_dz * Hp) + (Fp * dHp_dz) //Chain rule.

	Real const dFp_dz = ( 1.0 / ring_radius ) * (  - 2.0 * ( hp ) ) * ( 1 / ( sqrt ( rp3 ) ) );

	//EQ: dHp_dz = d_ellintk_dz + (dGp_dz * ellinte(mp) ) + ( Gp * d_ellinte_dz )
	//EQ: dHp_dz = (d_ellintk_dmp * dmp_dz )  + (dGp_dz * ellinte(mp) ) + ( Gp * dellinte_dmp * dmp_dz )

	Real const dmp_dz = ( 1.0 / ring_radius ) * ( ( rp * 0.0 ) - ( ( 4.0 * d ) * ( 2.0 * hp ) ) ) *  ( 1.0 / ( rp2 ) );  //Quotient Rule.

	Real const dGp_dz = ( 1.0 / ring_radius ) * ( ( sp * ( -2.0 * hp ) ) - ( tp * ( 2.0 * hp ) ) ) * ( 1.0 / ( sp2 ) );  //Quotient Rule.

	Real const dHp_dz = ( dellintk_dm( mp ) * dmp_dz )  + ( dGp_dz * ellinte( mp ) ) + ( Gp * dellinte_dm( mp ) * dmp_dz );

	dterm_dz = ( dFp_dz * Hp ) + ( Fp * dHp_dz );
}

///////////////////////////////////////////////////////////////

//atom_xyz,: the place where it should be calculated for
//molecular_ring_center : the middle of the plane where the atoms reside which carry the ring current
//base_z_axis: a vector perpendicular to the plane of the ring (perpendicular to the plane plane)

Real
delta_ring_current( numeric::xyzVector< core::Real > const & atom_xyz,
	numeric::xyzVector< core::Real > const & molecular_ring_center,
	numeric::xyzVector< core::Real > const & base_z_axis,
	RNA_CS_residue_parameters const & source_rsd_CS_params,
	Size const ring_ID )
{
	Real const rci = source_rsd_CS_params.ring_intensity( ring_ID );
	Real const rca = source_rsd_CS_params.ring_radius( ring_ID );
	Real const rch = source_rsd_CS_params.ring_height( ring_ID );

	Real rho = 0.0; //rho component of cylindrical coordinate
	Real z = 0.0;   //z   component of cylindrical coordinate

	get_rho_and_z( atom_xyz, molecular_ring_center, base_z_axis, rho, z );

	//Note definition of mterm and pterm is switched compared to NUCHEMIC!
	Real const mterm = delta_ring_current_term( rho, z, rca, (  - 1.0 ) * rch ); //This is contribution of the ring located at -rca below the molecular_plane
	Real const pterm = delta_ring_current_term( rho, z, rca, (  + 1.0 ) * rch ); //This is contribution of the ring located at +rca above the molecular_plane

	return ( -1.0 * source_rsd_CS_params.ring_current_coeff() * ( rci / rca ) * ( mterm + pterm ) ); //Definition uses opposite sign compare to NUCHEMICS!
}

///////////////////////////////////////////////////////////////
///The ring_current contribution of source_rsd to the chemical_shift at atom_xyz
Real
ring_current_effect( numeric::xyzVector< core::Real > const & atom_xyz, conformation::Residue const & source_rsd, RNA_CS_residue_parameters const & rna_cs_rsd_params ) {

	Real chem_shift = 0.0;

	if ( source_rsd.aa() != rna_cs_rsd_params.aa() ) utility_exit_with_message( "rsd.aa() != rna_cs_rsd_params.aa()!" );

	numeric::xyzMatrix< core::Real > const coordinate_matrix = get_rna_base_coordinate_system_from_CS_params( source_rsd, rna_cs_rsd_params );

	numeric::xyzVector< core::Real > const & base_z_axis = coordinate_matrix.col_z();

	for ( Size ring_ID = 1; ring_ID <= rna_cs_rsd_params.num_rings(); ring_ID++ ) {
		chem_shift += delta_ring_current( atom_xyz, ring_pos( source_rsd, rna_cs_rsd_params, ring_ID ), base_z_axis, rna_cs_rsd_params, ring_ID );
	}

	return chem_shift;
}


///////////////////////////////////////////////////////////////
//ONLY USE FOR TESTING PURPOSES!
Real
ring_current_effect_individual_ring( numeric::xyzVector< core::Real > const & atom_xyz, conformation::Residue const & source_rsd, RNA_CS_residue_parameters const & rna_cs_rsd_params, Size const source_ring_ID ){

	if ( source_rsd.aa() != rna_cs_rsd_params.aa() ) utility_exit_with_message( "rsd.aa() != rna_cs_rsd_params.aa()!" );

	numeric::xyzMatrix< core::Real > const coordinate_matrix = get_rna_base_coordinate_system_from_CS_params( source_rsd, rna_cs_rsd_params );

	numeric::xyzVector< core::Real > const & base_z_axis = coordinate_matrix.col_z();

	Real const chem_shift = delta_ring_current( atom_xyz, ring_pos( source_rsd, rna_cs_rsd_params, source_ring_ID ), base_z_axis, rna_cs_rsd_params, source_ring_ID );

	return chem_shift;
}
///////////////////////////////////////////////////////////////
//The gradient of ring_current_effect() with respect to r_vector (r_vector = CS_data_atom_xyz - molecular_ring_center)
//Between the source_ring_ID ring_current center of source_rsd and CS_data_atom_xyz.
//OK RIGHT NOW CONCERN WITH GETTING THIS FUNCTION RIGHT. OPTIMIZE LATER!

numeric::xyzVector< core::Real >
get_ring_current_deriv( numeric::xyzVector< core::Real > const & CS_data_atom_xyz,
	conformation::Residue const & source_rsd,
	core::Size const source_ring_ID,
	RNA_CS_residue_parameters const & source_rsd_CS_params ) {

	//Use cylindrical coordinates, ring_effect is a function of z and rho.
	//RC =( -1.0 * source_rsd_CS_params.ring_current_coeff() * (rci / rca) * (pterm(z, rho) + mterm(z, rho) ) );
	//Gradient  = dRC/drho * rho_norm + dRC/dz * z_norm

	Real const rci = source_rsd_CS_params.ring_intensity( source_ring_ID );
	Real const rca = source_rsd_CS_params.ring_radius( source_ring_ID );
	Real const rch = source_rsd_CS_params.ring_height( source_ring_ID );

	if ( source_rsd.aa() != source_rsd_CS_params.aa() ) utility_exit_with_message( "rsd.aa() != source_rsd_CS_params.aa()!" );

	numeric::xyzMatrix< core::Real > const coordinate_matrix = get_rna_base_coordinate_system_from_CS_params( source_rsd, source_rsd_CS_params );


	numeric::xyzVector< core::Real > const & base_z_axis = coordinate_matrix.col_z();

	numeric::xyzVector< core::Real > const & molecular_ring_center = ring_pos( source_rsd, source_rsd_CS_params, source_ring_ID );

	Real rho = 0.0; //rho component of cylindrical coordinate
	Real z = 0.0;   //z   component of cylindrical coordinate

	get_rho_and_z( CS_data_atom_xyz, molecular_ring_center, base_z_axis, rho, z );

	numeric::xyzVector< core::Real > const r_vector = CS_data_atom_xyz - molecular_ring_center;
	numeric::xyzVector< core::Real > const z_vector = ( base_z_axis * z );
	numeric::xyzVector< core::Real > const rho_vector = r_vector - ( base_z_axis * z );

	//Real const r_length = r_vector.length();

	numeric::xyzVector< core::Real > const z_norm = z_vector/z; //This gives back base_z_axis!;
	numeric::xyzVector< core::Real > const rho_norm = rho_vector/rho;

	/////////////////////////////////////////////////////////////////////

	Real dmterm_drho = 0.0;
	Real dmterm_dz  = 0.0;
	Real dpterm_drho = 0.0;
	Real dpterm_dz  = 0.0;

	get_ring_current_term_derivatives( rho, z, rca, (  - 1.0 ) * rch, dmterm_drho, dmterm_dz ); //This is contribution of the ring located at -rca below the molecular_plane

	get_ring_current_term_derivatives( rho, z, rca, (  + 1.0 ) * rch, dpterm_drho, dpterm_dz ); //This is contribution of the ring located at +rca above the molecular_plane

	Real const dRC_drho = -1.0 * source_rsd_CS_params.ring_current_coeff() * ( rci / rca ) * ( dmterm_drho + dpterm_drho );
	Real const dRC_dz   = -1.0 * source_rsd_CS_params.ring_current_coeff() * ( rci / rca ) * ( dmterm_dz + dpterm_dz );

	/////////////////////////////////////////////////////////////////////

	numeric::xyzVector< core::Real > const analytical_gradient = ( dRC_drho * rho_norm ) + ( dRC_dz * z_norm );

	return analytical_gradient;
}


} //chemical_shift
} //rna
} //scoring
} //core


