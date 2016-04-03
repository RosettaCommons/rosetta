// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file   core/scoring/rna/chemical_shift/RNA_CS_Util.cc
/// @brief
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#include <core/scoring/rna/chemical_shift/RNA_CS_Util.hh>

#include <ObjexxFCL/format.hh>
#include <core/chemical/AA.hh>

#include <numeric/xyzMatrix.hh>

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {


numeric::xyzMatrix< core::Real > const
get_rna_base_coordinate_system_from_CS_params( core::conformation::Residue const & rsd, RNA_CS_residue_parameters const & rna_cs_rsd_params ){ //x and y coordinate different from 'standard' get_rna_base_coordinate_system;

	using namespace chemical;

	if ( rsd.is_RNA() == false ) utility_exit_with_message( "rsd.is_RNA() == false!" );
	if ( rsd.aa() != rna_cs_rsd_params.aa() ) utility_exit_with_message( "rsd.aa() != rna_cs_rsd_params.aa()!" );

	numeric::xyzVector< core::Real > xf_xyz( 0.0, 0.0, 0.0 );
	numeric::xyzVector< core::Real > yf_xyz( 0.0, 0.0, 0.0 );
	numeric::xyzVector< core::Real > xt_xyz( 0.0, 0.0, 0.0 );
	numeric::xyzVector< core::Real > yt_xyz( 0.0, 0.0, 0.0 );

	Size xf_count = 0;
	Size yf_count = 0;

	Size xt_count = 0;
	Size yt_count = 0;

	Size const maxatoms = rna_cs_rsd_params.get_atomnames_size();

	for ( Size count = 1; count < maxatoms; count++ ) {

		Size const atom_index = rsd.atom_index( rna_cs_rsd_params.get_atomname( count ) );
		numeric::xyzVector< core::Real > const & atom_xyz = rsd.xyz( atom_index );

		if ( dround( rna_cs_rsd_params.atom_data( count, xdir ) ) == 1 ) {
			xf_count++;
			xf_xyz += atom_xyz;
		} else if ( dround( rna_cs_rsd_params.atom_data( count, xdir ) ) == 2 ) {
			xt_count++;
			xt_xyz += atom_xyz;
		}

		if ( dround( rna_cs_rsd_params.atom_data( count, ydir ) ) == 1 ) {
			yf_count++;
			yf_xyz += atom_xyz;
		} else if ( dround( rna_cs_rsd_params.atom_data( count, ydir ) ) == 2 ) {
			yt_count++;
			yt_xyz += atom_xyz;
		}
	}

	if ( xf_count != 1 && xf_count != 2 ) utility_exit_with_message( "xf_count != 1 && xf_count != 2" );
	if ( xt_count != 1 && xt_count != 2 ) utility_exit_with_message( "xt_count != 1 && xt_count != 2" );
	if ( yf_count != 1 && yf_count != 2 ) utility_exit_with_message( "yf_count != 1 && yf_count != 2" );
	if ( yt_count != 1 && yt_count != 2 ) utility_exit_with_message( "yt_count != 1 && yt_count != 2" );

	//Feb 21, 2012: UMM THIS FOLLOWS DEFINITION GIVEN IN NUCHEMICS..BUT x_base_coor and y_base_coor might not be completely perpendicular!

	numeric::xyzVector< core::Real > x_base_coor = ( ( xt_xyz / xt_count ) - ( xf_xyz / xf_count ) );
	x_base_coor.normalize();

	numeric::xyzVector< core::Real > y_base_coor = ( ( yt_xyz / yt_count ) - ( yf_xyz / yf_count ) );
	y_base_coor.normalize();

	numeric::xyzVector< core::Real > z_base_coor = cross( x_base_coor, y_base_coor );
	z_base_coor.normalize(); //necessary since x_base_coor and y_base_corr are not completely orthogonal!

	numeric::xyzMatrix< core::Real > const M = numeric::xyzMatrix< core::Real > ::cols( x_base_coor, y_base_coor, z_base_coor );
	return M;
}


/////////////////////////////////////////////////////////////////////////////


} //chemical_shift
} //rna
} //scoring
} //core
