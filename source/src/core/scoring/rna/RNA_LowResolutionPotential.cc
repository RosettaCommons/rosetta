// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_LowResolutionPotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>

// Package headers
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/interpolation_util.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>

#include <core/pose/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_ScoringUtil.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/io/izstream.hh>

#include <ObjexxFCL/FArray2A.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/interpolation/periodic_range/full/interpolation.hh>

// C++

#include <basic/Tracer.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


static basic::Tracer tr( "core.scoring.rna.RNA_LowResolutionPotential" );

namespace core {
namespace scoring {
namespace rna {

core::Real RNA_LowResolutionPotential::dummy_deriv = 0.0;

typedef  numeric::xyzMatrix< Real > Matrix;

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

// Hey should we define a copy constructor and a clone()?

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
RNA_LowResolutionPotential::RNA_LowResolutionPotential():
	////////////////////////////////////////////
	rna_basepair_radius_cutoff_( 8.0 ),
	rna_basepair_stagger_cutoff_( 3.0 ),
	rna_basepair_radius_cutoff2_( rna_basepair_radius_cutoff_ * rna_basepair_radius_cutoff_ ), // 64.0
	basepair_xy_bin_width_( 2.0 ),
	basepair_xy_num_bins_( 10 ),
	basepair_xy_table_size_( 10 ),
	basepair_xy_z_fade_zone_( 0.5 ),
	////////////////////////////////////////////
	base_stack_min_height_( 2.4 ),
	base_stack_max_height_( 6.0 ),
	base_stack_radius_( 4.0 ),
	base_stack_radius2_( base_stack_radius_ * base_stack_radius_ ),
	base_stack_z_fade_zone_( 0.5 ),
	base_stack_rho_fade_zone_( 0.5 ),
	////////////////////////////////////////////
	axis_bin_width_( 0.2 ),
	axis_num_bins_ ( 11  ),
	////////////////////////////////////////////
	stagger_num_bins_( 11 ),
	stagger_bin_width_( 0.4 ),
	stagger_distance_cutoff_( 2.0 ),
	////////////////////////////////////////////
	base_backbone_bin_width_( 1.0 ),
	base_backbone_num_bins_( 16 ),
	base_backbone_table_size_( 8 ),
	base_backbone_distance_cutoff_ ( 12.0 ),
	base_backbone_z_cutoff_ ( 2.0 ),
	base_backbone_rho_cutoff_( 8.0 ),
	base_backbone_atom_dist_cutoff_( 4.0 ),
	base_backbone_z_fade_zone_( 0.25 ),
	base_backbone_rho_fade_zone_( 0.5 ),
	base_backbone_check_atom_neighbor_( false ), /*doesn't seem compatible with deriv*/
	////////////////////////////////////////////
	backbone_backbone_bin_width_( 0.25 ),
	backbone_backbone_distance_cutoff_ ( 6.0 ),
	backbone_backbone_num_bins_ ( 20 ),
	////////////////////////////////////////////
	rna_repulsive_max_penalty_( 8.0 ),
	rna_repulsive_screen_scale_( 2.5 ),
	rna_repulsive_distance_cutoff_( 8.0 ),
	rna_repulse_all_( true ),
	num_RNA_base_pair_orientations_( 2 ), /*parallel/anti*/
	num_RNA_backbone_oxygen_atoms_( 6 ), /*backbone oxygens*/
	num_RNA_res_types_( 4 ), /*a/c/g/u*/
	o2prime_index_within_special_backbone_atoms_( 6 ),
	o2p_index_within_special_backbone_atoms_( 2 ),
	interpolate_( true ), //Turn this off to match Rosetta++, up to bug fixes; turn it on to allow correct derivative calculation.
	fade_( true ), //Turn this off to match Rosetta++; needed to prevent hard boundaries in base pairing + stacking potentials.
	rna_verbose_( false ),
	more_precise_base_pair_classification_( false )
{
	// These don't *have* to be hard-wired numbers.
	// if we're more clever about the file read-in.
	// For now, just copy what was in Rosetta++.
	rna_basepair_xy_.dimension( basepair_xy_num_bins_, basepair_xy_num_bins_, num_RNA_base_pair_orientations_, num_RNA_res_types_, 4 /*a/c/g/u*/ );
	rna_axis_.dimension( axis_num_bins_ );
	rna_stagger_.dimension( stagger_num_bins_ );
	rna_base_backbone_xy_.dimension( base_backbone_num_bins_, base_backbone_num_bins_,
																	 num_RNA_res_types_, num_RNA_backbone_oxygen_atoms_ );
	rna_backbone_backbone_potential_.dimension( backbone_backbone_num_bins_ );
	rna_repulsive_weight_.dimension( num_RNA_backbone_oxygen_atoms_ );

	assert( Size( basepair_xy_bin_width_   * basepair_xy_num_bins_/ 2.0 )   == basepair_xy_table_size_ );
	assert( Size( base_backbone_bin_width_ * base_backbone_num_bins_/ 2.0 ) == base_backbone_table_size_ );

	rna_basepair_xy_ = 0.0;
	rna_axis_ = 0.0;
	rna_stagger_ = 0.0;
	rna_base_backbone_xy_ = 0.0;
	rna_backbone_backbone_potential_ = 0.0;
	rna_repulsive_weight_ = 0.0;

	initialize_rna_basepair_xy();
	initialize_rna_axis();
	initialize_rna_stagger();
	initialize_RNA_backbone_oxygen_atoms();
	initialize_atom_numbers_for_backbone_score_calculations();
	initialize_rna_base_backbone_xy();
	initialize_rna_backbone_backbone_weights();
	initialize_rna_backbone_backbone();
	initialize_rna_repulsive_weights();
	initialize_more_precise_base_pair_cutoffs();

}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_rna_basepair_xy(){

	std::string const filename( "scoring/rna/rna_base_pair_xy.dat" );
  utility::io::izstream data_stream(  basic::database::full_name( filename )  );

	if ( !data_stream ) {
		std::cerr << "Can't find specified basepair potential file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
		return;
	}

	tr << "Reading basepair x - y potential file: " << filename << std::endl; ;
	// read data
	Size res1, res2, xbin, ybin, direction;
	Real potential;
	while ( data_stream >> xbin ) {
		data_stream >> ybin >> res1 >> res2 >> direction >> potential >> skip ;
		rna_basepair_xy_( xbin, ybin, res1, res2, direction ) = potential;
	}

	tr << "Finished reading basepair x - y potential file: " << filename << std::endl; ;

	//close file
	data_stream.close();
}

////////////////////////////////////////////////////////////////////////////////////////
// Uh, should put this in a file.
void
RNA_LowResolutionPotential::initialize_rna_axis(){
   rna_axis_(  1 ) = -3.319;
   rna_axis_(  2 ) = -0.938;
   rna_axis_(  3 ) = 0.000;
   rna_axis_(  4 ) = 0.000;
   rna_axis_(  5 ) = 0.000;
   rna_axis_(  6 ) = -0.309;
   rna_axis_(  7 ) = -0.748;
   rna_axis_(  8 ) = -1.559;
   rna_axis_(  9 ) = -2.576;
   rna_axis_( 10 ) = -3.478;
   rna_axis_( 11 ) = -3.529;
}


///////////////////////////////////////////////////////////////////////////////////////
// Uh, should put this in a file.
void
RNA_LowResolutionPotential::initialize_rna_stagger(){
	rna_stagger_(  1 ) = 0.000;
	rna_stagger_(  2 ) = 0.000;
	rna_stagger_(  3 ) = -0.062;
	rna_stagger_(  4 ) = -1.174;
	rna_stagger_(  5 ) = -2.085;
	rna_stagger_(  6 ) = -2.897;
	rna_stagger_(  7 ) = -2.496;
	rna_stagger_(  8 ) = -2.253;
	rna_stagger_(  9 ) = -1.156;
	rna_stagger_( 10 ) = -0.461;
	//	rna_stagger_( 11 ) = -0.354;
	// This needs to asymptote to zero, or derivatives don't look so nice.
	rna_stagger_( 11 ) = 0.000;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_RNA_backbone_oxygen_atoms(){

	RNA_backbone_oxygen_atoms_.clear();
	RNA_backbone_oxygen_atoms_.push_back( " OP2" );
	RNA_backbone_oxygen_atoms_.push_back( " OP1" );
	RNA_backbone_oxygen_atoms_.push_back( " O5'" );
	RNA_backbone_oxygen_atoms_.push_back( " O4'" );
	RNA_backbone_oxygen_atoms_.push_back( " O3'" );
	RNA_backbone_oxygen_atoms_.push_back( " O2'" );

	//Useful for preventing string lookups:
	o2p_index_within_special_backbone_atoms_    =  2;
	o2prime_index_within_special_backbone_atoms_ =  6;


}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_rna_base_backbone_xy(){

	std::string const filename ( "scoring/rna/rna_base_backbone_xy.dat" );

	// open file
  utility::io::izstream data_stream( basic::database::full_name( filename ) );

	if ( !data_stream ) {
		std::cerr << "Can't find specified non - base - base potential file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
		return;
	}

	tr << "Reading non - base - base x - y potential file: " << filename << std::endl; ;
	// read data
	Size res1, xbin, ybin, atomindex;
	Real potential;
	while ( data_stream >> xbin ) {
		data_stream >> ybin >> res1 >> atomindex >> potential >> skip ;
		rna_base_backbone_xy_( xbin, ybin, res1, atomindex ) = potential;
	}

	//close file
	data_stream.close();
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_rna_backbone_backbone_weights(){
	rna_backbone_backbone_weight_.dimension( 6 );
	rna_backbone_backbone_weight_ = 0.0;

	rna_backbone_backbone_weight_( 1 ) = 1.0; // OP2
	rna_backbone_backbone_weight_( 2 ) = 1.0; // OP1
	rna_backbone_backbone_weight_( 6 ) = 2.0; // O2'
}


//////////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_rna_backbone_backbone()
{

	std::string const filename = "scoring/rna/rna_backbone_backbone.dat";

	// open file
  utility::io::izstream data_stream(  basic::database::full_name( filename )  );

	if ( !data_stream ) {
		std::cerr << "Can't find specified RNA backbone - backbone potential file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
		return;
	}

	tr << "Reading RNA backbone backbone potential file: " << filename << std::endl; ;
	// read data
	Size rbin;
	Real r1, r2, potential;
	while ( data_stream >> r1 ) {
		data_stream >> r2 >> rbin >> potential >> skip ;
		rna_backbone_backbone_potential_( rbin ) = potential;
	}

	//close file
	data_stream.close();
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_rna_repulsive_weights(){

	// Note that unless specified by user "-rna_phosphate_repulse_all"
	//  only OP1-OP1 repulsion will be calculated. See eval_rna_phosphate_score().
	rna_repulsive_weight_( 1 ) = 1.0; // OP2
	rna_repulsive_weight_( 2 ) = 1.0; // OP1
	rna_repulsive_weight_( 3 ) = 1.0; // O5'
	rna_repulsive_weight_( 4 ) = 1.0; // O4'
	rna_repulsive_weight_( 5 ) = 1.0; // O3'
	// O2' --> weight stays at zero.
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::get_rna_basepair_xy(
    Distance const x,
		Distance const y,
		Distance const z, // used in fading
		Real const cos_theta,
		conformation::Residue const & res_i,
		conformation::Residue const & res_j,
		bool const update_deriv /* = true */,
		Real & deriv_x,
		Real & deriv_y,
		Real & deriv_z ) const
{
	Size const theta_bin = ( cos_theta < 0 ) ?  1 : 2;

	assert( res_i.is_RNA() );
	assert( res_j.is_RNA() );

	Size const res_i_bin = core::chemical::rna::convert_acgu_to_1234( res_i.name1() );
	Size const res_j_bin = core::chemical::rna::convert_acgu_to_1234( res_j.name1() );

	assert( res_i_bin > 0 );
	assert( res_j_bin > 0 );

	Real value( 0.0 );
	deriv_x = 0.0;
	deriv_y = 0.0;
	deriv_z = 0.0;

	if ( interpolate_ ) {
		Real interpolate_value( 0.0 );
		FArray2A < Real > ::IR const zero_index( 0, basepair_xy_num_bins_ - 1 );
		ObjexxFCL::FArray2A < Real > const rna_basepair_xy_for_res_res( rna_basepair_xy_( 1, 1, theta_bin, res_i_bin, res_j_bin ), zero_index, zero_index );
		using namespace numeric::interpolation::periodic_range::full;

		if ( update_deriv ) {
			interpolate_value = bilinearly_interpolated(
				 Real( x + basepair_xy_table_size_ ),
				 Real( y + basepair_xy_table_size_ ),
				 basepair_xy_bin_width_,
				 basepair_xy_num_bins_,
				 rna_basepair_xy_for_res_res,
				 deriv_x, deriv_y );
		} else {
			interpolate_value = bilinearly_interpolated(
				 Real( x + basepair_xy_table_size_ ),
				 Real( y + basepair_xy_table_size_ ),
				 basepair_xy_bin_width_,
				 basepair_xy_num_bins_,
				 rna_basepair_xy_for_res_res );
		}
		//		std::cout << "BASE_BASE: " << x << " " << y << " ==> " << value <<  " " << interpolate_value << std::endl;

		value = interpolate_value;

	} else { // old Rosetta++-style, no interpolation
		Size x_bin, y_bin;

		x_bin = Size( ( x + basepair_xy_table_size_ ) / basepair_xy_bin_width_ ) + 1;
		y_bin = Size( ( y + basepair_xy_table_size_ ) / basepair_xy_bin_width_ ) + 1;

		if ( x_bin < 1 ) x_bin = 1;
		if ( y_bin < 1 ) y_bin = 1;
		if ( x_bin > basepair_xy_table_size_ ) x_bin = basepair_xy_table_size_;
		if ( y_bin > basepair_xy_table_size_ ) y_bin = basepair_xy_table_size_;

		value =  rna_basepair_xy_( x_bin, y_bin, theta_bin, res_i_bin, res_j_bin );

	}

	if ( fade_ ) {
		Real fade_value( 1.0 ), fade_deriv( 0.0 );

		//First do z correction
		get_fade_correction( z, -1.0 * rna_basepair_stagger_cutoff_, rna_basepair_stagger_cutoff_,
												      basepair_xy_z_fade_zone_ /* 0.5 */, fade_value, fade_deriv );

		//		std::cout << "BASE_BASE: " << x << " " << y << " ==> " << value <<  " " << fade_value;

		deriv_x  = deriv_x * fade_value;
		deriv_y  = deriv_y * fade_value;
		deriv_z  = value * fade_deriv;
		value *= fade_value;

		// Now rho fading ( radial ).
 		Distance const rho = std::sqrt( x*x + y*y );
 		if ( rho > 0.0 ) { //better be!
			get_fade_correction( rho, -rna_basepair_radius_cutoff_, rna_basepair_radius_cutoff_,
 													 base_stack_rho_fade_zone_ /* 0.5 */, fade_value, fade_deriv );

			//			std::cout << " RHO_FADE  " << fade_value << std::endl;

 			deriv_x = fade_value * deriv_x  +  value * fade_deriv * ( x / rho ) ;
			deriv_y = fade_value * deriv_y  +  value * fade_deriv * ( y / rho );
 			deriv_z *= fade_value;
 			value   *= fade_value;

 		}

	}

	return value;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::get_rna_stack_score(
		Distance const x, // used in fading
		Distance const y, // used in fading
		Distance const z, // used in fading
		Real & deriv_x,
		Real & deriv_y,
		Real & deriv_z ) const
{
	// In Rosetta++, this was just a constant.
	Real value( -0.5 );
	deriv_x = 0.0;
	deriv_y = 0.0;
	deriv_z = 0.0;

	if ( fade_ ) {
		// First do z-fading.
		Real fade_value( 1.0 ), fade_deriv( 0.0 );
		if ( z > 0 ) {
			get_fade_correction( z, base_stack_min_height_, base_stack_max_height_,
													 base_stack_z_fade_zone_ /* 0.5 */, fade_value, fade_deriv );
		} else {
			get_fade_correction( z, -1 * base_stack_max_height_, -1 * base_stack_min_height_,
													 base_stack_z_fade_zone_ /* 0.5 */, fade_value, fade_deriv );
		}

		deriv_z = value * fade_deriv;
		value   *= fade_value;

		// Now rho fading ( radial ).
 		Distance const rho = std::sqrt( x*x + y*y );
 		if ( rho > 0.0 ) { //better be!
			get_fade_correction( rho, -base_stack_radius_, base_stack_radius_,
 													 base_stack_rho_fade_zone_ /* 0.5 */, fade_value, fade_deriv );
 			deriv_x = value * fade_deriv * ( x / rho ) ;
			deriv_y = value * fade_deriv * ( y / rho );
 			deriv_z *= fade_value;
 			value   *= fade_value;
 		}

	}

	return value;
}



///////////////////////////////////////////////////////////////////////////////
//Following could be part of its own class...
Real
RNA_LowResolutionPotential::get_rna_axis_score(
   Real const cos_theta,
	 Real & deriv ) const{

	Size cos_theta_bin;
	// cos_theta_bin = static_cast<int> ( ( cos_theta + 1.0f)/ axis_bin_width_ ) + 1;
	// Note that above, which was in Rosetta++ is slightly off.
	// Bin boundaries are [-1.1,-0.9], [-0.9,-0.7], ... [0.9,1.1]. So should be:
	cos_theta_bin = static_cast< int > ( ( cos_theta + 1.1f )/ axis_bin_width_ ) + 1;

	Real value ( 0.0 );

	if ( interpolate_ ) {
		Real interpolate_value( 0.0 );
		//interpolate_value_and_deriv( rna_axis_, axis_bin_width_, cos_theta + 1.0, interpolate_value, deriv );
		interpolate_value_and_deriv( rna_axis_, axis_bin_width_, cos_theta + 1.1, interpolate_value, deriv );
		value = interpolate_value;
		//std::cout << "AXIS: " << F(8,3,cos_theta) << " " << F(8,3,value) << "  vs. " << F(8,3,interpolate_value ) << std::endl;
	} else {
		assert( cos_theta_bin > 0 && cos_theta_bin <= axis_num_bins_ );
		value = rna_axis_( cos_theta_bin );
	}
	return value;
}

///////////////////////////////////////////////////////////////////////////////
//Following could be part of its own class...
Real
RNA_LowResolutionPotential::get_rna_stagger_score( Distance const height, Real & deriv ) const
{

	Real value( 0.0 );
	deriv = 0.0;


	//////////////////////////////////////////////////////////////////////////////////////////
	// SILLY HACK -- no the problem isn't here.
	//deriv = 2 * height;
	//return height * height;
	//////////////////////////////////////////////////////////////////////////////////////////

	if ( std::abs( height ) > stagger_distance_cutoff_ ) return 0.0;

	if ( interpolate_ ) {
			interpolate_value_and_deriv( rna_stagger_, stagger_bin_width_, height + 2.2, value, deriv );
			//	std::cout << "STAGGER: " << F(8,3,height) << " " << F(8,3,value) << " " << F(8,3,interpolate_value) << std::endl;
	} else {
		//Following, from Rosetta++, is slightly off. Bins are at [-2.2,-1.8],[-1.8,-1.6],...[1.8,2.2]
		//	Size height_bin = static_cast<int> ( (height + 2.0)/ stagger_bin_width_ ) + 1;
		Size height_bin = static_cast< int > ( ( height + 2.2 )/ stagger_bin_width_ ) + 1;
		if ( height_bin < 1 || height_bin > stagger_num_bins_ ) height_bin = 1;
		value = rna_stagger_( height_bin );
	}

	return value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::get_rna_base_backbone_xy(
	Distance const x,
	Distance const y,
	Distance const z,
	conformation::Residue const & res_i,
	Size const & atom_num_j_bin,
	bool const update_deriv /* = false */,
    Real & deriv_x /* = dummy_deriv */,
    Real & deriv_y /* = dummy_deriv */,
	Real & deriv_z /* = dummy_deriv */ ) const
{

	Real value( 0.0 );
	deriv_x = 0.0;
	deriv_y = 0.0;
	deriv_z = 0.0;

	Size const res_i_bin = core::chemical::rna::convert_acgu_to_1234( res_i.name1() );

	assert( res_i_bin > 0 );
	assert( atom_num_j_bin > 0 );

	if ( interpolate_ ) {
		//		Real interpolate_value( 0.0 );
		FArray2A < Real > ::IR const zero_index( 0, base_backbone_num_bins_ - 1 );
		ObjexxFCL::FArray2A < Real > const rna_base_backbone_xy_for_res_atom( rna_base_backbone_xy_( 1, 1, res_i_bin, atom_num_j_bin ), zero_index, zero_index );
		using namespace numeric::interpolation::periodic_range::full;

		if ( update_deriv || true ) {
			value = bilinearly_interpolated(
				Real( x + base_backbone_table_size_ ),
				Real( y + base_backbone_table_size_ ),
				base_backbone_bin_width_,
				base_backbone_num_bins_,
				rna_base_backbone_xy_for_res_atom,
				deriv_x, deriv_y );
		} else {
			value = bilinearly_interpolated(
				Real( x + base_backbone_table_size_ ),
				Real( y + base_backbone_table_size_ ),
				base_backbone_bin_width_,
				base_backbone_num_bins_,
				rna_base_backbone_xy_for_res_atom );
		}
		//		value = interpolate_value;
	} else { //old school.
		Size x_bin, y_bin;
		x_bin = Size( ( x + base_backbone_table_size_ ) ) + 1;
		y_bin = Size( ( y + base_backbone_table_size_ ) ) + 1;

		if ( x_bin < 1 ) x_bin = 1;
		if ( y_bin < 1 ) y_bin = 1;
		if ( x_bin > 2 * base_backbone_table_size_ ) x_bin = 2 * base_backbone_table_size_;
		if ( y_bin > 2 * base_backbone_table_size_ ) y_bin = 2 * base_backbone_table_size_;

		value = rna_base_backbone_xy_( x_bin, y_bin, res_i_bin, atom_num_j_bin );
	}


	if ( fade_ ) {
		Real fade_value( 1.0 ), fade_deriv( 0.0 );

		//First do z correction
		get_fade_correction( z, -1.0 * base_backbone_z_cutoff_, base_backbone_z_cutoff_,
												      base_backbone_z_fade_zone_ /* 0.5 */, fade_value, fade_deriv );

		deriv_x  = deriv_x * fade_value;
		deriv_y  = deriv_y * fade_value;
		deriv_z  = value * fade_deriv;
		value *= fade_value;

		//		std::cout << "BASE_BACKBONE_Z_FADE: " << z << " " << base_backbone_z_cutoff_  <<  " " << base_backbone_z_fade_zone_ << " ==> " << fade_value << std::endl;

 		// Now rho fading ( radial ).
		Distance const rho = std::sqrt( x*x + y*y );
		if ( rho > 0.0 ) { //better be!
 			get_fade_correction( rho, -1.0 * base_backbone_rho_cutoff_, base_backbone_rho_cutoff_,
													 base_backbone_rho_fade_zone_ /* 0.5 */, fade_value, fade_deriv );

			deriv_x = fade_value * deriv_x  +  value * fade_deriv * ( x / rho ) ;
 			deriv_y = fade_value * deriv_y  +  value * fade_deriv * ( y / rho );
			deriv_z *= fade_value;
			value   *= fade_value;

			//			std::cout << "BASE_BACKBONE_RHO_FADE: " << rho << " " << base_backbone_rho_cutoff_ <<  " " << base_backbone_rho_fade_zone_ << " ==> " << fade_value << std::endl;
 		}

	}


	//	std::cout << "BASE_BACKBONE: " << x << " " << y << " ==> " << value <<  " " << interpolate_value << std::endl;

	return value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::get_rna_backbone_backbone_score(
	Distance const & r,
	Size const & atom_num_j_bin,
	Real & deriv /* = 0.0 */
) const
{

	Real value( 0.0 );
	deriv = 0.0;

	if ( interpolate_ ) {
		interpolate_value_and_deriv( rna_backbone_backbone_potential_, backbone_backbone_bin_width_, r, value, deriv );
	} else {
		Size const bin = static_cast< Size > ( r/ backbone_backbone_bin_width_ ) + 1;
		if ( bin <= backbone_backbone_num_bins_ )		value = rna_backbone_backbone_potential_( bin );
	}

	value *= rna_backbone_backbone_weight_( atom_num_j_bin );
	deriv *= rna_backbone_backbone_weight_( atom_num_j_bin );

	//	if ( value < 0.0 ) {
	//		std::cout  << "RNA_BACKBONE_BACKBONE: " << value << ".  R: " << r  << "  ATOMNUM:  " <<
	//			atom_num_j_bin << std::endl;
	//	}

	return value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::get_rna_repulsive_score(
   Distance const & r,
	 Size const & atom_num_j_bin,
	 Real & deriv /* = 0.0 */ ) const
{

	static Real const rna_repulsive_max_penalty_ = 8.0; //In kT.
	static Real const rna_repulsive_screen_scale_ = 2.5; //In Angstroms
	static Real const rna_repulsive_distance_cutoff_ = 8.0; //In Angstroms

	static Real offset = rna_repulsive_max_penalty_ * exp( -1.0 * rna_repulsive_distance_cutoff_ / rna_repulsive_screen_scale_ );

	deriv = 0.0;
	if ( r < rna_repulsive_distance_cutoff_ ){
		Real potential = rna_repulsive_max_penalty_ * exp( -1.0 * r/ rna_repulsive_screen_scale_ ) - offset;
		deriv = -1.0 * ( rna_repulsive_max_penalty_ / rna_repulsive_screen_scale_ ) *
			exp( -1.0 * r/ rna_repulsive_screen_scale_ );

		potential *=  rna_repulsive_weight_( atom_num_j_bin );
		deriv     *=  rna_repulsive_weight_( atom_num_j_bin );

		return potential;
	}

	return 0.0;

}


///////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::update_rna_centroid_info(
	pose::Pose & pose
) const
{
	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//Doesn't recalculate stuff if already updated:
	rna_centroid_info.update( pose );
}

///////////////////////////////////////////////////////////////////////////////
// DEPRECATED?
void
RNA_LowResolutionPotential::update_rna_base_base_interactions(
	pose::Pose & pose
) const
{

	// If we want to play some cache-ing tricks later...
	//	FArray2D_bool const & pair_moved( pose.get_pair_moved() );

	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//Doesn't recalculate stuff if already updated:
	rna_centroid_info.update( pose );

	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );


	//	Real const Z_CUTOFF( 2.5 );

	Size const total_residue = pose.total_residue();

	rna::RNA_RawBaseBaseInfo & rna_raw_base_base_info( rna_scoring_info.rna_raw_base_base_info() );
	rna_raw_base_base_info.resize( total_residue );
	//	if (rna_raw_base_base_info.calculated()) return;

	ObjexxFCL::FArray3D < Real > & base_pair_array( rna_raw_base_base_info.base_pair_array() );
	ObjexxFCL::FArray3D < Real > & base_axis_array( rna_raw_base_base_info.base_axis_array() );
	ObjexxFCL::FArray3D < Real > & base_stagger_array( rna_raw_base_base_info.base_stagger_array() );
	ObjexxFCL::FArray2D < Real > & base_stack_array( rna_raw_base_base_info.base_stack_array() );
	ObjexxFCL::FArray2D < Real > & base_stack_axis_array( rna_raw_base_base_info.base_stack_axis_array() );
	ObjexxFCL::FArray2D < Real > & base_geometry_orientation_array( rna_raw_base_base_info.base_geometry_orientation_array() );
	ObjexxFCL::FArray2D < Real > & base_geometry_height_array( rna_raw_base_base_info.base_geometry_height_array() );

	//Following may change in the future, if we
	// decide to take advantage of "pair_moved" arrays to
	// prevent recalculation of energy terms.
	base_pair_array = 0.0;
	base_stagger_array = 0.0;
	base_axis_array = 0.0;
	base_stack_array = 0.0;
	base_stack_axis_array = 0.0;
	base_geometry_orientation_array = 0.0;
	base_geometry_height_array = 0.0;

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	//Main loop
	// NOTE: Following evaluates each base pair from both sides, because base pair
	// terms are asymmetric in i<->j. So loop is over all neighbors!
	for ( Size i = 1; i <= total_residue; i++ ){
		conformation::Residue const & res_i( pose.residue( i ) );

		if ( !res_i.is_RNA() ) continue;

		//Note that this centroid and stubs could be calculated once at the beginning of the scoring!!!
		Vector const & centroid_i( base_centroids[i] );
		kinematics::Stub const & stub_i( base_stubs[i] );
    Matrix const & M_i( stub_i.M );
    Vector const & x_i = M_i.col_x();
    Vector const & y_i = M_i.col_y();
    Vector const & z_i = M_i.col_z();

		for ( graph::Graph::EdgeListConstIter
					 iter = energy_graph.get_node( i )->const_edge_list_begin();
				 iter != energy_graph.get_node( i )->const_edge_list_end();
				 ++iter ){

			Size j( ( *iter )->get_other_ind( i ) );

			//			if ( !pair_moved(i,j) && rna_array_state_ok ) continue;

			conformation::Residue const & res_j( pose.residue( j ) );
			if ( !res_j.is_RNA() ) continue;

			Vector const & centroid_j( base_centroids[j] );
			kinematics::Stub const & stub_j( base_stubs[j] );

			Vector d_ij = centroid_j - centroid_i;
			Real const dist_x = dot_product( d_ij, x_i );
			Real const dist_y = dot_product( d_ij, y_i );
			Real const dist_z = dot_product( d_ij, z_i );
			Real const rho2 = dist_x*dist_x + dist_y*dist_y;

			//This could all be precomputed, actually.
			Matrix const & M_j( stub_j.M );
			//			Vector const & x_j = M_j.col_x();
			//			Vector const & y_j = M_j.col_y();
			Vector const & z_j = M_j.col_z();
			Real const cos_theta = dot_product( z_i, z_j );
			//This arc-cosine costs a lot actually:
			//			Real const theta = numeric::conversions::degrees( numeric::arccos( cos_theta ) );

			//Is it a base-pair, a base-stack, or do we ignore it?
			Size edge_bin( 1 );
			if ( ( std::abs( dist_z ) < rna_basepair_stagger_cutoff_ ) ){
				//A possible base pair
				if ( rho2 < rna_basepair_radius_cutoff2_ ){
					//BASE PAIR
					Real temp_rna_bp_score( 0.0 );

					temp_rna_bp_score = get_rna_basepair_xy( dist_x, dist_y, dist_z, cos_theta, res_i, res_j );

					// Some of the base pairing scores are 0.0. By default, pad by
					// a tiny tiny amount to make sure that they're still
					// counted as very weak base pairs during later book-keeping.
					temp_rna_bp_score -= 0.0001;

					Real zeta_hoogsteen_cutoff( 60.0 ), zeta_sugar_cutoff( -60.0 );
					get_zeta_cutoff( res_i, zeta_hoogsteen_cutoff, zeta_sugar_cutoff );

					Real const zeta = numeric::conversions::degrees( std::atan2( dist_y, dist_x ) );
					if ( zeta < zeta_hoogsteen_cutoff && zeta > zeta_sugar_cutoff )      edge_bin = core::chemical::rna::WATSON_CRICK;  //Watson-Crick edge
					else if ( zeta > zeta_hoogsteen_cutoff )   edge_bin = core::chemical::rna::HOOGSTEEN; // Hoogsteen edge
					else                       edge_bin = core::chemical::rna::SUGAR; // Sugar edge

					if ( rna_verbose_ ){
						Real const theta = numeric::conversions::degrees( numeric::arccos( cos_theta ) );
						tr        << " Possible base pair: "
											<< res_i.name3() << I( 3, i ) << "-"
											<< res_j.name3() << I( 3, j )
											<< " edge: " << I( 1, edge_bin )
											<< "  dists " << F( 4, 2, dist_x )
											<< " " << F( 4, 2, dist_y )
											<< " " << F( 4, 2, dist_z )
											<< " theta " << F( 5, 1, theta )
											<< "; rho " << F( 4, 2, std::sqrt( rho2 ) )
											<< "; zeta " << F( 4, 2, zeta )
											<< "; zeta_cut " << F( 4, 2, zeta_hoogsteen_cutoff )
											<< " : SCORE " << F( 6, 4, temp_rna_bp_score )
							//												<< " " << get_rna_axis_score( theta) << " " << get_rna_stagger_score( dist_z )
											<<	std::endl;
					}

					base_pair_array    ( i, j, edge_bin ) = temp_rna_bp_score;
					base_axis_array    ( i, j, edge_bin ) = get_rna_axis_score( cos_theta );
					base_stagger_array ( i, j, edge_bin ) = get_rna_stagger_score( dist_z );

				}
			}

			if ( std::abs( dist_z ) >= base_stack_min_height_ &&
					 std::abs( dist_z ) <= base_stack_max_height_  &&
					 rho2 < base_stack_radius2_ ) {
				//Possible BASE STACK1
				base_stack_array ( i, j ) = get_rna_stack_score( dist_x, dist_y, dist_z );
				base_stack_axis_array ( i, j ) = get_rna_axis_score( cos_theta );
			} //basepair

			base_geometry_orientation_array ( i, j ) = cos_theta;
			base_geometry_height_array      ( i, j ) = dist_z;

		} //j
	} //i

	//	rna_raw_base_base_info.set_calculated( true );

}

void
RNA_LowResolutionPotential::update_rna_base_pair_list(
	pose::Pose & pose
) const
{
	using namespace core::scoring::rna;

	pose.update_residue_neighbors();
	update_rna_base_base_interactions( pose );

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_RawBaseBaseInfo & raw_base_base_info( rna_scoring_info.rna_raw_base_base_info() );
	rna::RNA_FilteredBaseBaseInfo & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	// Maybe we should put an if statement.
	rna_filtered_base_base_info.carry_out_filtering( raw_base_base_info );
}

///////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::eval_rna_base_pair_energy(
  rna::RNA_RawBaseBaseInfo & rna_raw_base_base_info,
  conformation::Residue const & rsd1,
  conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Vector const & centroid1,
	Vector const & centroid2,
	kinematics::Stub const & stub1,
	kinematics::Stub const & stub2
) const
{
	// Following fills in arrays inside rna_raw_base_base_info, which you have to
	// look inside to get computed energies.
	eval_rna_base_pair_energy_one_way( rna_raw_base_base_info, rsd1, rsd2, pose, centroid1, centroid2, stub1, stub2 );
	eval_rna_base_pair_energy_one_way( rna_raw_base_base_info, rsd2, rsd1, pose, centroid2, centroid1, stub2, stub1 );
}

/////////////////////////////////////////////////////////////
// These are slightly more sane cutoffs for defining what is
// the Watson-Crick edge vs. Hoogsteen edge vs. sugar edge.
void
RNA_LowResolutionPotential::setup_precise_zeta_cutoffs( chemical::AA const & na_rad,
																												std::string const & hoogsteen_cutoff_atom,
																												std::string const & sugar_cutoff_atom
)
{
	using namespace core::chemical;
	using namespace core::conformation;

	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_RNA );
	RNA_CentroidInfo rna_centroid_info;

	ResidueOP rsd = ResidueFactory::create_residue( *( rsd_set->aa_map( na_rad ) )[1] );
	Vector const & centroid_i = rna_centroid_info.get_base_centroid( *rsd );
	kinematics::Stub const & stub_i = rna_centroid_info.get_base_coordinate_system( *rsd, centroid_i );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();

	Size const res_i_bin = core::chemical::rna::convert_acgu_to_1234( rsd->name1() );

	Real dist_x( 0.0 ), dist_y( 0.0 ), zeta( 0.0 );

	Vector d_ij_hoogsteen = rsd->xyz( hoogsteen_cutoff_atom ) - centroid_i;
	dist_x = dot_product( d_ij_hoogsteen, x_i );
	dist_y = dot_product( d_ij_hoogsteen, y_i );
	zeta = numeric::conversions::degrees( std::atan2( dist_y, dist_x ) );
	zeta_hoogsteen_cutoff_precise_( res_i_bin ) = zeta;

	if ( true ) {
		Vector d_ij_sugar = rsd->xyz( sugar_cutoff_atom ) - centroid_i;
		dist_x = dot_product( d_ij_sugar, x_i );
		dist_y = dot_product( d_ij_sugar, y_i );
		zeta = numeric::conversions::degrees( std::atan2( dist_y, dist_x ) );
		zeta_sugar_cutoff_precise_( res_i_bin ) = zeta + 10.0;
	}

	//zeta_hoogsteen_cutoff_precise_( res_i_bin ) = 	zeta_hoogsteen_cutoff_precise_( res_i_bin ) - 120.0;

	//std::cout << "HEY! PRECISION CUTOFF FOR ZETA: " <<
	//		zeta_hoogsteen_cutoff_precise_( res_i_bin ) << " " <<
	//		zeta_sugar_cutoff_precise_( res_i_bin ) << " " <<  std::endl;

}

////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_more_precise_base_pair_cutoffs() //Doesn't do anything if already initialized.
{

	using namespace core::chemical;

	zeta_hoogsteen_cutoff_precise_.dimension ( 4 );
	zeta_sugar_cutoff_precise_.dimension ( 4 );

	setup_precise_zeta_cutoffs( na_rad, " H61", " C2 " );
	setup_precise_zeta_cutoffs( na_rcy, " H42", " O2 " );
	setup_precise_zeta_cutoffs( na_rgu, " O6 ", " N2 " );
	setup_precise_zeta_cutoffs( na_ura, " O4 ", " O2 " );


}

////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::get_zeta_cutoff(
  conformation::Residue const & res_i,
	Real & zeta_hoogsteen_cutoff,
	Real & zeta_sugar_cutoff
) const
{
	if ( more_precise_base_pair_classification_ ) {
		Size const res_i_bin = core::chemical::rna::convert_acgu_to_1234( res_i.name1() );
		zeta_hoogsteen_cutoff = zeta_hoogsteen_cutoff_precise_( res_i_bin );
		zeta_sugar_cutoff = zeta_sugar_cutoff_precise_( res_i_bin );
	} else {
		zeta_hoogsteen_cutoff = 60.0;
		zeta_sugar_cutoff = -60.0;
	}
}

///////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::eval_rna_base_pair_energy_one_way(
  rna::RNA_RawBaseBaseInfo & rna_raw_base_base_info,
  conformation::Residue const & res_i,
  conformation::Residue const & res_j,
	pose::Pose const & pose,
	Vector const & centroid_i,
	Vector const & centroid_j,
	kinematics::Stub const & stub_i,
	kinematics::Stub const & stub_j
) const
{

	if ( !res_i.is_RNA() ) return;
	if ( !res_j.is_RNA() ) return;

	Size const total_residue = pose.total_residue();

	// Is this necessary?
	rna_raw_base_base_info.resize( total_residue );

	ObjexxFCL::FArray3D < Real > & base_pair_array( rna_raw_base_base_info.base_pair_array() );
	ObjexxFCL::FArray3D < Real > & base_axis_array( rna_raw_base_base_info.base_axis_array() );
	ObjexxFCL::FArray3D < Real > & base_stagger_array( rna_raw_base_base_info.base_stagger_array() );
	ObjexxFCL::FArray2D < Real > & base_stack_array( rna_raw_base_base_info.base_stack_array() );
	ObjexxFCL::FArray2D < Real > & base_stack_axis_array( rna_raw_base_base_info.base_stack_axis_array() );
	ObjexxFCL::FArray2D < Real > & base_geometry_orientation_array( rna_raw_base_base_info.base_geometry_orientation_array() );
	ObjexxFCL::FArray2D < Real > & base_geometry_height_array( rna_raw_base_base_info.base_geometry_height_array() );

	Size const i( res_i.seqpos() );
	Size const j( res_j.seqpos() );

	// Zero out these arrays for the residue pair of interest.
	for ( Size k = 1; k <= core::chemical::rna::NUM_EDGES; k++ ){
	  base_pair_array( i, j, k ) = 0.0;
	  base_axis_array( i, j, k ) = 0.0;
	  base_stagger_array( i, j, k ) = 0.0;
	}
	base_stack_array( i, j ) = 0.0;
	base_stack_axis_array( i, j ) = 0.0;
	base_geometry_orientation_array( i, j ) = 0.0;
	base_geometry_height_array( i, j ) = 0.0;

	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	Vector d_ij = centroid_j - centroid_i;
	Real const dist_x = dot_product( d_ij, x_i );
	Real const dist_y = dot_product( d_ij, y_i );
	Real const dist_z = dot_product( d_ij, z_i );
	Real const rho2 = dist_x*dist_x + dist_y*dist_y;

	//This could all be precomputed, actually.
	Matrix const & M_j( stub_j.M );
	//			Vector const & x_j = M_j.col_x();
	//			Vector const & y_j = M_j.col_y();
	Vector const & z_j = M_j.col_z();
	Real const cos_theta = dot_product( z_i, z_j );
	//This arc-cosine costs a lot actually:
	//			Real const theta = numeric::conversions::degrees( numeric::arccos( cos_theta ) );

	//Is it a base-pair, a base-stack, or do we ignore it?
	Size edge_bin( 1 );
	if ( ( std::abs( dist_z ) < rna_basepair_stagger_cutoff_ ) ){
		//A possible base pair
		if ( rho2 < rna_basepair_radius_cutoff2_ ){
			//BASE PAIR
			Real temp_rna_bp_score( 0.0 );

			temp_rna_bp_score = get_rna_basepair_xy( dist_x, dist_y, dist_z, cos_theta, res_i, res_j );

			// Some of the base pairing scores are 0.0. By default, pad by
			// a tiny tiny amount to make sure that they're still
			// counted as very weak base pairs during later book-keeping.
			temp_rna_bp_score -= 0.0001;

			Real zeta_hoogsteen_cutoff( 60.0 ), zeta_sugar_cutoff( -60.0 );
			get_zeta_cutoff( res_i, zeta_hoogsteen_cutoff, zeta_sugar_cutoff );

			Real const zeta = numeric::conversions::degrees( std::atan2( dist_y, dist_x ) );
			if ( zeta < zeta_hoogsteen_cutoff && zeta > zeta_sugar_cutoff )      edge_bin = core::chemical::rna::WATSON_CRICK;  //Watson-Crick edge
			else if ( zeta > zeta_hoogsteen_cutoff )   edge_bin = core::chemical::rna::HOOGSTEEN; // Hoogsteen edge
			else                       edge_bin = core::chemical::rna::SUGAR; // Sugar edge

			if ( rna_verbose_ ){
				Real const theta = numeric::conversions::degrees( numeric::arccos( cos_theta ) );
				tr        << " Possible base pair: "
									<< res_i.name3() << I( 3, i ) << "-"
									<< res_j.name3() << I( 3, j )
									<< " edge: " << I( 1, edge_bin )
									<< "  dists " << F( 4, 2, dist_x )
									<< " " << F( 4, 2, dist_y )
									<< " " << F( 4, 2, dist_z )
									<< " Centroid " << centroid_i( 1 ) << " " << centroid_i( 2 ) << " " << centroid_i( 3 )
									<< " CENTROID " << centroid_j( 1 ) << " " << centroid_j( 2 ) << " " << centroid_j( 3 )
									<< " theta " << F( 5, 1, theta )
									<< "; rho " << F( 4, 2, std::sqrt( rho2 ) )
									<< " : SCORE " << F( 6, 4, temp_rna_bp_score )
					//												<< " " << get_rna_axis_score( theta) << " " << get_rna_stagger_score( dist_z )
									<<	std::endl;
			}

			base_pair_array    ( i, j, edge_bin ) = temp_rna_bp_score;
			base_axis_array    ( i, j, edge_bin ) = get_rna_axis_score( cos_theta );
			base_stagger_array ( i, j, edge_bin ) = get_rna_stagger_score( dist_z );

		}
	}

	if ( std::abs( dist_z ) >= base_stack_min_height_ &&
			 std::abs( dist_z ) <= base_stack_max_height_  &&
			 rho2 < base_stack_radius2_ ) {
		//Possible BASE STACK1
		base_stack_array ( i, j ) = get_rna_stack_score( dist_x, dist_y, dist_z );
		base_stack_axis_array ( i, j ) = get_rna_axis_score( cos_theta );
	} //basepair

	base_geometry_orientation_array ( i, j ) = cos_theta;
	base_geometry_height_array      ( i, j ) = dist_z;

}


///////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::eval_atom_derivative_base_base(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
		scoring::EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2 ) const
{

	Size const & i = atom_id.rsd();
	conformation::Residue const & res_i( pose.residue( i ) );
	Size const & atom_num_i( atom_id.atomno() );

	if ( !res_i.is_RNA() ) return;

	//First an easy filter -- only need to put derivs on base's torsion.

	if ( atom_num_i != core::chemical::rna::chi1_torsion_atom_index( res_i ) ) return;

	// Information saved from the last score.
	rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	//	rna::RNA_RawBaseBaseInfo const & raw_base_base_info( rna_scoring_info.rna_raw_base_base_info() );
	rna::RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );

	//	assert( rna_filtered_base_base_info.calculated()  );

	// Axis and stagger terms will have hard boundaries unless they're
	// multiplied by base pair (x,y) or base stack (z) terms.
	assert( rna_filtered_base_base_info.scale_axis_stagger() );

	// Yea, this is what we need...
  ObjexxFCL::FArray2D < Real > const & filtered_base_pair_array ( rna_filtered_base_base_info.filtered_base_pair_array() );
  ObjexxFCL::FArray2D < Real > const & filtered_base_axis_array ( rna_filtered_base_base_info.filtered_base_axis_array() );
  ObjexxFCL::FArray2D < Real > const & filtered_base_stagger_array( rna_filtered_base_base_info.filtered_base_stagger_array() );
  ObjexxFCL::FArray2D < Real > const & filtered_base_stack_array( rna_filtered_base_base_info.filtered_base_stack_array() );
  ObjexxFCL::FArray2D < Real > const & filtered_base_stack_axis_array( rna_filtered_base_base_info.filtered_base_stack_axis_array() );

	rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//	assert( rna_centroid_info.calculated()  );
	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	//	Vector const & heavy_atom_i_xyz( res_i.xyz( atom_num_i ) ) ;

	Vector const & centroid_i( base_centroids[i] );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( graph::Graph::EdgeListConstIter
				 iter = energy_graph.get_node( i )->const_edge_list_begin();
			 iter != energy_graph.get_node( i )->const_edge_list_end();
			 ++iter ){

		Size j( ( *iter )->get_other_ind( i ) );

		conformation::Residue const & res_j( pose.residue( j ) );

		Vector const & centroid_j( base_centroids[j] );
		kinematics::Stub const & stub_j( base_stubs[j] );

		Matrix const & M_j( stub_j.M );
		Vector const & x_j = M_j.col_x();
		Vector const & y_j = M_j.col_y();
		Vector const & z_j = M_j.col_z();
		Real const cos_theta = dot_product( z_i, z_j );
		//This arc-cosine costs a lot actually:
		//		Real const theta = numeric::conversions::degrees( numeric::arccos( cos_theta ) );

		//First cycle through all base pairs that count... anything involving this residue?
		if ( filtered_base_pair_array( i, j ) < 0.0 ) {

			Real const base_axis_score_unscaled    = filtered_base_axis_array( i, j )/filtered_base_pair_array( i, j );
			Real const base_stagger_score_unscaled = filtered_base_stagger_array( i, j )/filtered_base_pair_array( i, j );
			Real bp_deriv_x( 0.0 ), bp_deriv_y( 0.0 ), bp_deriv_z( 0.0 ), axis_deriv( 0.0 ), stagger_deriv( 0.0 );
			Real const & basepair_axis_stagger_scaling( rna_filtered_base_base_info.basepair_axis_stagger_scaling() );

			// If its been scored, all the various geometry cutoffs should already have been satisfied
			{ // from this base coordinate system to other base.
				Vector d_ij = centroid_j - centroid_i;
				Real const dist_x = dot_product( d_ij, x_i );
				Real const dist_y = dot_product( d_ij, y_i );
				Real const dist_z = dot_product( d_ij, z_i );
				//		Real const rho2 = dist_x*dist_x + dist_y*dist_y;

				get_rna_basepair_xy( dist_x, dist_y, dist_z, cos_theta, res_i, res_j,
														true /*update_deriv*/, bp_deriv_x, bp_deriv_y, bp_deriv_z );
				get_rna_axis_score( cos_theta, axis_deriv );
				get_rna_stagger_score( dist_z, stagger_deriv );

				//Time for the chain rule... to get the total force and torque.
				// 1. Correct for scale_axis_stagger craziness.
				// 2. I don't fully understand factor of two yet.
				Vector const f2 =
					( x_i * bp_deriv_x + y_i * bp_deriv_y + z_i * bp_deriv_z )
					* ( weights[ rna_base_pair ] +
							weights[ rna_base_axis ] * base_axis_score_unscaled +
							weights[ rna_base_stagger ] *	base_stagger_score_unscaled )
					- z_i * basepair_axis_stagger_scaling * filtered_base_pair_array( i, j ) * weights[ rna_base_stagger ] * stagger_deriv;

				Vector const f1 = cross( f2, centroid_j )
					- basepair_axis_stagger_scaling * filtered_base_pair_array( i, j ) * weights[ rna_base_axis ] * axis_deriv * cross( z_i, z_j );

				F1 -= f1;
				F2 -= f2;

				//				std::cout << "BASEPAIR_DERIV1: " << i << " " << j << " " << cos_theta << " " << f1(1) << " " << f2(1) << std::endl;
			}

			{ // from other base coordinate system to this one.
				Vector d_ji = centroid_i - centroid_j;
				Real const dist_x = dot_product( d_ji, x_j );
				Real const dist_y = dot_product( d_ji, y_j );
				Real const dist_z = dot_product( d_ji, z_j );
				//		Real const rho2 = dist_x*dist_x + dist_y*dist_y;

				get_rna_basepair_xy( dist_x, dist_y, dist_z, cos_theta, res_j, res_i,
														true /*update_deriv*/, bp_deriv_x, bp_deriv_y, bp_deriv_z );
				get_rna_axis_score( cos_theta, axis_deriv );
				get_rna_stagger_score( dist_z, stagger_deriv );

				//Time for the chain rule... to get the total force and torque.
				// Correct for scale_axis_stagger craziness.
				Vector const f2 =
					( x_j * bp_deriv_x + y_j * bp_deriv_y + z_j * bp_deriv_z )
					* ( weights[ rna_base_pair ] +
							weights[ rna_base_axis ] * base_axis_score_unscaled +
							weights[ rna_base_stagger ] *	base_stagger_score_unscaled )
					- z_j * basepair_axis_stagger_scaling * filtered_base_pair_array( j, i ) * weights[ rna_base_stagger ] * stagger_deriv;

				Vector const f1 = cross( f2, centroid_i )
					- basepair_axis_stagger_scaling * filtered_base_pair_array( j, i ) * weights[ rna_base_axis ] * axis_deriv * cross( z_j, z_i );

				F1 += f1;
				F2 += f2;

				//				std::cout << "BASEPAIR_DERIV2: " << i << " " << j << " " << cos_theta << " " << f1(1) << " " << f2(1) << std::endl;
			}

		} else if ( filtered_base_stack_array( i, j ) < 0.0 ) {
			//Then cycle through all base stacks that count.

			Real const base_stack_axis_score_unscaled    = filtered_base_stack_axis_array( i, j )/filtered_base_stack_array( i, j );
			Real bs_deriv_x( 0.0 ), bs_deriv_y( 0.0 ), bs_deriv_z( 0.0 ), axis_deriv( 0.0 );
			Real const & basestack_axis_scaling( rna_filtered_base_base_info.basestack_axis_scaling() );

			{
				Vector d_ij = centroid_j - centroid_i;
				Real const dist_x = dot_product( d_ij, x_i );
				Real const dist_y = dot_product( d_ij, y_i );
				Real const dist_z = dot_product( d_ij, z_i );

				// If its been scored, all the various geometry cutoffs should already have been satisfied!
				get_rna_stack_score( dist_x, dist_y, dist_z, bs_deriv_x, bs_deriv_y, bs_deriv_z );
				get_rna_axis_score( cos_theta, axis_deriv );

				Vector const f2 =
					( x_i * bs_deriv_x + y_i * bs_deriv_y + z_i * bs_deriv_z ) *
					( weights[ rna_base_stack ] +
						weights[ rna_base_stack_axis ] * base_stack_axis_score_unscaled );
				Vector const f1 = cross( f2, centroid_j )
					- basestack_axis_scaling * weights[ rna_base_stack_axis ] * filtered_base_stack_array( i, j ) * axis_deriv * cross( z_i, z_j );

				F1 -= f1;
				F2 -= f2;

				//				std::cout << "BASESTACK_DERIV1: " << i << " " << j << " " << cos_theta << " " << f1(1) << " " << f2(1) << std::endl;
			}
			{
				Vector d_ji = centroid_i - centroid_j;
				Real const dist_x = dot_product( d_ji, x_j );
				Real const dist_y = dot_product( d_ji, y_j );
				Real const dist_z = dot_product( d_ji, z_j );

				// If its been scored, all the various geometry cutoffs should already have been satisfied!
				get_rna_stack_score( dist_x, dist_y, dist_z, bs_deriv_x, bs_deriv_y, bs_deriv_z );
				get_rna_axis_score( cos_theta, axis_deriv );

				//				std::cout << "HELLOSTACK --> " << dist_x << " " << dist_y << " " << dist_z << " " << bs_deriv_x << " " << bs_deriv_y << " " << bs_deriv_y << " " << weights[ rna_base_stack_axis ] << " " << base_stack_axis_score_unscaled << std::endl;

				Vector const f2 =
					( x_j * bs_deriv_x + y_j * bs_deriv_y + z_j * bs_deriv_z ) *
					( weights[ rna_base_stack ] +
						weights[ rna_base_stack_axis ] * base_stack_axis_score_unscaled );
				Vector const f1 = cross( f2, centroid_i )
					- basestack_axis_scaling * weights[ rna_base_stack_axis ] * filtered_base_stack_array( j, i ) * axis_deriv * cross( z_j, z_i );

				F1 += f1;
				F2 += f2;

				//				std::cout << "BASESTACK_DERIV2: " << i << " " << j << " " << cos_theta << " " << f1(1) << " " << f2(1) << std::endl;
			}
		}

	}

}


//////////////////////////////////////////////////////////////////////////////////
// This actually imposes a hard-edge, and in special cases this causes a
// problem with (numerical) derivative calculations, and accumulation at the boundary.
// So I added a weight that is very very close to one unless we're on the edge (within 0.01 A),
// in which case we fade...
// In principle we should also put in the derivative at the edge... but not yet implemented.
bool
RNA_LowResolutionPotential::check_for_base_neighbor(
  conformation::Residue const & rsd1,
	Vector const & heavy_atom_j,
	Real & atom_cutoff_weight ) const
{

	atom_cutoff_weight = 1.0;
	if ( !base_backbone_check_atom_neighbor_ ) return true;

	Size const num_heavy_atoms ( rsd1.nheavyatoms() );
	Size const rna_base_start  ( rsd1.first_sidechain_atom() );

	bool found_a_neighbor = false;
	Real min_nbr_dist2( 1000.0 );
	for ( Size atom_num_i = rna_base_start; atom_num_i <= num_heavy_atoms; atom_num_i++ ){
		Vector const & base_heavy_atom_i( rsd1.xyz( atom_num_i ) );
		Real const nbr_dist2 = ( heavy_atom_j - base_heavy_atom_i ).length_squared();
		if ( nbr_dist2 < min_nbr_dist2 ){
			min_nbr_dist2 = nbr_dist2;
			//			break;
		}
	}

	Real const min_nbr_dist = std::sqrt( min_nbr_dist2 );
	if ( min_nbr_dist < base_backbone_atom_dist_cutoff_ ) found_a_neighbor = true;

	Real deriv_currently_ignored( 0.0 );
	if ( found_a_neighbor ) {
		get_fade_correction( min_nbr_dist,
												 -base_backbone_atom_dist_cutoff_ /*will not come into play*/,
												 base_backbone_atom_dist_cutoff_,
												 0.02 /*fade boundary*/,
												 atom_cutoff_weight,
												 deriv_currently_ignored );
		//		std::cout << "Hello? NBR_DIST" << nbr_dist << " " << atom_cutoff_weight << std::endl;
	}

	return found_a_neighbor;
}

/////////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::rna_base_backbone_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Vector const & centroid1,
	Vector const & centroid2,
	kinematics::Stub const & stub1,
	kinematics::Stub const & stub2
) const
{
	return ( rna_base_backbone_pair_energy_one_way( rsd1, rsd2, centroid1, stub1 ) +
					 rna_base_backbone_pair_energy_one_way( rsd2, rsd1, centroid2, stub2  ) );
}

//////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::rna_base_backbone_pair_energy_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Vector const & centroid_i,
	kinematics::Stub const & stub_i
) const
{
	////////////////////////////////////////////////////////////////////////////////
	//Hmm, is rsd1 always less than rsd2?

	if ( !rsd1.is_RNA() ) return 0.0;
	if ( !rsd2.is_RNA() ) return 0.0;

	if ( rsd1.is_coarse() ) return 0.0;  //coarse-grained!
	if ( rsd2.is_coarse() ) return 0.0;  //coarse-grained!

	Size const i = rsd1.seqpos();
	Size const j = rsd2.seqpos();

	if ( abs( static_cast< int > ( i - j ) ) < 2 ) return 0.0;

	//	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	//	rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//	assert( rna_centroid_info.calculated() ); //This needs to all be calculated earlier, dude.

	//	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
	//  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );
	//	Vector const & centroid_i( base_centroids[i] );
	//	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	Real total_score( 0.0 );

	//For speed, a cached set of atom numbers that go with RNA_backbone_oxygen_atoms_ (which is a bunch of strings)
	//	ObjexxFCL::FArray2D< Size > const &
	//		atom_numbers_for_backbone_score_calculations = rna_scoring_info.atom_numbers_for_backbone_score_calculations();

	// Go over sugar and phosphate oxygen atoms
	for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ){

		//		std::string const atom_j = RNA_backbone_oxygen_atoms_[ m ];
		Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ]; //atom_numbers_for_backbone_score_calculations( rsd2.seqpos(),  m );

		Vector const heavy_atom_j( rsd2.xyz( atom_num_j ) );

		Vector const d_ij = heavy_atom_j - centroid_i;

		Real const dist_ij = d_ij.length();

		if ( dist_ij < base_backbone_distance_cutoff_ ){

			Real const dist_x = dot_product( d_ij, x_i );
			Real const dist_y = dot_product( d_ij, y_i );
			Real const dist_z = dot_product( d_ij, z_i );

			Real const rho = std::sqrt( dist_x * dist_x + dist_y * dist_y );

			if ( std::abs( dist_z ) > base_backbone_z_cutoff_ ) continue; // Look for atoms in the base plane
			if ( rho > base_backbone_rho_cutoff_ ) continue; // Look for atoms in the base plane

			//sanity check...
			//make sure we're in H-bonding distance of some base atom.
			Real atom_cutoff_weight( 1.0 );
			if ( !check_for_base_neighbor( rsd1, heavy_atom_j, atom_cutoff_weight ) ) continue;

			Real const score_contribution =
				atom_cutoff_weight * get_rna_base_backbone_xy( dist_x, dist_y, dist_z, rsd1, m );

		  total_score += score_contribution;

			if ( rna_verbose_ ){
				tr <<
					"BASE - BACKBONE " <<
					rsd1.name3() <<
					I( 3, i ) << " " <<
					rsd2.atom_name( atom_num_j )  <<  " " <<
					I( 3, j ) << " " <<
					" [" << F( 4, 2, rho ) << ", " << F( 4, 2, dist_z ) << "]:  " <<
							F( 6, 2, score_contribution ) <<
					std::endl;
			}

		}

	}

	return total_score;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lot of code copying here -- should minirosetta++ derivative calculation be refactored?
//   Not so much for speed but due to code redundancy in defining loops, etc.
//
void
RNA_LowResolutionPotential::eval_atom_derivative_rna_base_backbone(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		Vector & F1,
 		Vector & F2 ) const
{

	Size const & i = atom_id.rsd();
	conformation::Residue const & rsd1( pose.residue( i ) );
	Size const & atom_num_i( atom_id.atomno() );

	F1 = 0.0;
	F2 = 0.0;

	if ( !rsd1.is_RNA() ) return;
	assert( interpolate_ );

	//For speed, a cached set of atom numbers that go with RNA_backbone_oxygen_atoms_ (which is a bunch of strings)
	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	//	ObjexxFCL::FArray2D< Size > const &
	//		atom_numbers_for_backbone_score_calculations = rna_scoring_info.atom_numbers_for_backbone_score_calculations();
	rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );

	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	//Either we're a backbone oxygen atom, or we're the RNA first base atom (N1 or N9).

	//	std::cout << rsd1.aa() << " CHI1 TORSION ATOM ==> "  << chi1_torsion_atom_index( rsd1 )  << " " << chi1_torsion_atom( rsd1 )  << std::endl;
	if ( atom_num_i == core::chemical::rna::chi1_torsion_atom_index( rsd1 ) ) {

		Vector const & centroid_i( base_centroids[i] );
		kinematics::Stub const & stub_i( base_stubs[i] );
		Matrix const & M_i( stub_i.M );
		Vector const & x_i = M_i.col_x();
		Vector const & y_i = M_i.col_y();
		Vector const & z_i = M_i.col_z();

		//		Size const num_heavy_atoms ( rsd1.nheavyatoms() );
		//		Size const rna_base_start( rsd1.first_sidechain_atom() );

		for ( graph::Graph::EdgeListConstIter
				iter = energy_graph.get_node( i )->const_edge_list_begin();
				iter != energy_graph.get_node( i )->const_edge_list_end();
				++iter ){

			Size j( ( *iter )->get_other_ind( i ) );

			if ( abs( static_cast< int > ( i - j ) ) < 2 ) continue;

			conformation::Residue const & rsd2( pose.residue( j ) );
			if ( !rsd2.is_RNA() ) continue;

			// Go over sugar and phosphate oxygen atoms
			for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ){

				//		std::string const atom_j = RNA_backbone_oxygen_atoms_[ m ];
				Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ]; //atom_numbers_for_backbone_score_calculations( rsd2.seqpos(),  m );

				Vector const heavy_atom_j( rsd2.xyz( atom_num_j ) );

				Vector const d_ij = heavy_atom_j - centroid_i;

				Real const dist_ij = d_ij.length();

				if ( dist_ij >= base_backbone_distance_cutoff_ ) continue;

				Real const dist_x = dot_product( d_ij, x_i );
				Real const dist_y = dot_product( d_ij, y_i );
				Real const dist_z = dot_product( d_ij, z_i );

				Real const rho = std::sqrt( dist_x * dist_x + dist_y * dist_y );

				if ( std::abs( dist_z ) > base_backbone_z_cutoff_ ) continue; // Look for atoms in the base plane
				if ( rho > base_backbone_rho_cutoff_ ) continue; // Look for atoms in the base plane

				//sanity check...
				//make sure we're in H-bonding distance of some base atom.
				Real atom_cutoff_weight( 1.0 );
				if ( !check_for_base_neighbor( rsd1, heavy_atom_j, atom_cutoff_weight ) ) continue;

				Real deriv_x( 0.0 ), deriv_y( 0.0 ), deriv_z( 0.0 );
				get_rna_base_backbone_xy( dist_x, dist_y, dist_z, rsd1, m, true /*update_deriv*/, deriv_x, deriv_y, deriv_z );

				Vector const f2 = -1.0 * ( x_i * deriv_x + y_i * deriv_y + z_i * deriv_z );
				Vector const f1 = cross( f2, heavy_atom_j );

				F1 += f1;
				F2 += f2;

				//				std::cout << "BACKBONEBASE_DERIV1 " << f2(1) << std::endl;

			} // m
		} // nbrs
	} else {

		Size n = find_backbone_oxygen_atom( atom_num_i );

		Vector const heavy_atom_i( rsd1.xyz( atom_num_i ) );

		if ( n > 0 ){

			for ( graph::Graph::EdgeListConstIter
						 iter = energy_graph.get_node( i )->const_edge_list_begin();
					 iter != energy_graph.get_node( i )->const_edge_list_end();
					 ++iter ){

				Size j( ( *iter )->get_other_ind( i ) );

				if ( abs( static_cast< int > ( i - j ) ) < 2 ) continue;

				conformation::Residue const & rsd2( pose.residue( j ) );
				if ( !rsd2.is_RNA() ) continue;

				Vector const & centroid_j( base_centroids[j] );
				kinematics::Stub const & stub_j( base_stubs[j] );
				Matrix const & M_j( stub_j.M );
				Vector const & x_j = M_j.col_x();
				Vector const & y_j = M_j.col_y();
				Vector const & z_j = M_j.col_z();

				Vector const d_ij = heavy_atom_i - centroid_j;

				Real const dist_ij = d_ij.length();

				if ( dist_ij >= base_backbone_distance_cutoff_ ) continue;

				Real const dist_x = dot_product( d_ij, x_j );
				Real const dist_y = dot_product( d_ij, y_j );
				Real const dist_z = dot_product( d_ij, z_j );

				Real const rho = std::sqrt( dist_x * dist_x + dist_y * dist_y );

				if ( std::abs( dist_z ) > base_backbone_z_cutoff_ ) continue; // Look for atoms in the base plane
				if ( rho > base_backbone_rho_cutoff_ ) continue; // Look for atoms in the base plane

				//sanity check...
				//make sure we're in H-bonding distance of some base atom.
				Real atom_cutoff_weight( 1.0 );
				if ( !check_for_base_neighbor( rsd2, heavy_atom_i, atom_cutoff_weight ) ) continue;

				Real deriv_x( 0.0 ), deriv_y( 0.0 ), deriv_z( 0.0 );

				get_rna_base_backbone_xy( dist_x, dist_y, dist_z, rsd2, n, true /*update_deriv*/, deriv_x, deriv_y, deriv_z );

				Vector const f2 = x_j * deriv_x + y_j * deriv_y + z_j * deriv_z;
				Vector const f1 = cross( f2, heavy_atom_i );

				F1 += f1;
				F2 += f2;

				//				std::cout << "BACKBONEBASE_DERIV2 " << f1(1) << std::endl;

			} // n
		} // nbrs

	}

	//	std::cout << "BACKBONEBASE_DERIV_SUM " << F1(1) << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::rna_backbone_backbone_pair_energy(
			conformation::Residue const & rsd1,
			conformation::Residue const & rsd2
) const
{
	return ( rna_backbone_backbone_pair_energy_one_way( rsd1, rsd2 ) +
					 rna_backbone_backbone_pair_energy_one_way( rsd2, rsd1 ) ) ;
}

///////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::rna_backbone_backbone_pair_energy_one_way(
			conformation::Residue const & rsd1,
			conformation::Residue const & rsd2
) const
{

	if ( !rsd1.is_RNA() ) return 0.0;
	if ( !rsd2.is_RNA() ) return 0.0;

	if ( rsd1.is_coarse() ) return 0.0;  //coarse-grained!
	if ( rsd2.is_coarse() ) return 0.0;  //coarse-grained!

	Real rna_backbone_backbone_score = 0.0;

	//  static Real const dist_cutoff ( 6.0 );

	//	std::string const atom_i = " O2'";

	Size const i( rsd1.seqpos() );
	Size const j( rsd2.seqpos() );

	if ( abs( static_cast< int > ( i - j ) ) <= 2 ) return 0.0;

	//For speed, a cached set of atom numbers that go with RNA_backbone_oxygen_atoms_ (which is a bunch of strings)
	//	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	//	ObjexxFCL::FArray2D< Size > const &
	//		atom_numbers_for_backbone_score_calculations = rna_scoring_info.atom_numbers_for_backbone_score_calculations();

	Size const atom_num_i = atom_numbers_for_backbone_score_calculations_[ o2prime_index_within_special_backbone_atoms_ ];
	Vector const & heavy_atom_i = rsd1.xyz( atom_num_i );

	// Go over sugar and phosphate oxygen atoms!
	for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ){

		if ( rna_backbone_backbone_weight_( m ) < 0.001 ) continue;

		Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ];

		//Don't double-count 2'-OH <--> 2'-OH interactions
		if ( atom_num_j == atom_num_i && j < i ) continue;

		Vector const & heavy_atom_j = rsd2.xyz( atom_num_j );
		Vector const d_ij = heavy_atom_j - heavy_atom_i;
		Distance const dist_ij = d_ij.length();

		if ( dist_ij > backbone_backbone_distance_cutoff_ ) continue;

		Real const score_contribution = get_rna_backbone_backbone_score( dist_ij, m );
		rna_backbone_backbone_score += score_contribution;

		if ( rna_verbose_ && score_contribution < -0.01 ){
			std::string const atom_j = RNA_backbone_oxygen_atoms_[ m ];
			tr <<
				"BACKBONE_BACKBONE " <<
				rsd1.name1() <<
				I( 3, i ) << " " <<
				rsd1.atom_name( atom_num_i ) <<
				" " << rsd1.xyz( atom_num_i )[1] <<
				" -- " <<
				rsd2.name1() <<
				I( 3, j ) << " " <<
				rsd2.atom_name( atom_num_j ) << " " <<
				" " << rsd2.xyz( atom_num_j )[1] <<
				F( 6, 2, score_contribution ) << " [" << dist_ij << "]" <<
						std::endl;
		}

	}

	return rna_backbone_backbone_score;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_LowResolutionPotential::find_backbone_oxygen_atom(
   ObjexxFCL::FArray2D < Size > const & atom_numbers_for_backbone_score_calculations,
	 Size const & i,
	 Size const & atom_num_i ) const
{

	Size n( 0 );
	bool is_backbone_oxygen_atom( false );

	for ( n = 1; n <= num_RNA_backbone_oxygen_atoms_; n++ ){
		if ( atom_numbers_for_backbone_score_calculations( i, n ) == atom_num_i ) {
			is_backbone_oxygen_atom = true; break;
		}
	}
	if ( !is_backbone_oxygen_atom ) return 0;

	return n;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_LowResolutionPotential::find_backbone_oxygen_atom(
	 Size const & atom_num_i ) const
{

	Size n( 0 );
	bool is_backbone_oxygen_atom( false );

	for ( n = 1; n <= num_RNA_backbone_oxygen_atoms_; n++ ){
		if ( atom_numbers_for_backbone_score_calculations_[ n ] == atom_num_i ) {
			is_backbone_oxygen_atom = true; break;
		}
	}
	if ( !is_backbone_oxygen_atom ) return 0;

	return n;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lot of code copying here -- should minirosetta++ derivative calculation be refactored?
//   Not so much for speed but due to code redundancy in defining loops, etc.
//
void
RNA_LowResolutionPotential::eval_atom_derivative_rna_backbone_backbone(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		Vector & F1,
 		Vector & F2 ) const
{

	Size const & i = atom_id.rsd();
	conformation::Residue const & rsd1( pose.residue( i ) );
	Size const & atom_num_i( atom_id.atomno() );

	F1 = 0.0;
	F2 = 0.0;

	if ( !rsd1.is_RNA() ) return;

	assert( interpolate_ );

	//For speed, a cached set of atom numbers that go with RNA_backbone_oxygen_atoms_ (which is a bunch of strings)
	//	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	//	ObjexxFCL::FArray2D< Size > const &
	//		atom_numbers_for_backbone_score_calculations = rna_scoring_info.atom_numbers_for_backbone_score_calculations();

	// FIRST WAY, check if this atom is OP1, cycle over other oxygen atoms.
	Size const atom_num_o2prime = atom_numbers_for_backbone_score_calculations_[  o2prime_index_within_special_backbone_atoms_ ];

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	if ( atom_num_i == atom_num_o2prime ) {
		// FIRST WAY, check if this atom is 2'-OH, cycle over other oxygen atoms.

		Vector const & heavy_atom_i = rsd1.xyz( atom_num_i );

		for ( graph::Graph::EdgeListConstIter
				iter = energy_graph.get_node( i )->const_edge_list_begin();
				iter != energy_graph.get_node( i )->const_edge_list_end();
				++iter ){

			Size j( ( *iter )->get_other_ind( i ) );

			if ( abs( static_cast< int > ( i - j ) ) <= 2 ) continue;

			conformation::Residue const & rsd2( pose.residue( j ) );
			if ( !rsd2.is_RNA() ) continue;

			// Go over sugar and phosphate oxygen atoms!
			for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ){

				if ( rna_backbone_backbone_weight_( m ) < 0.001 ) continue;

				Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ];

				//Don't double-count 2'-OH <--> 2'-OH interactions
				//			if ( atom_num_j == atom_num_i && j < i) continue;

				Vector const & heavy_atom_j = rsd2.xyz( atom_num_j );
				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Distance const dist_ij = d_ij.length();

				if ( dist_ij > backbone_backbone_distance_cutoff_ ) continue;

				Real deriv( 0.0 );
				get_rna_backbone_backbone_score( dist_ij, m, deriv );

				Vector const f2( -1.0 * deriv * d_ij/ dist_ij );
				Vector const f1(  1.0 * cross( f2, heavy_atom_j ) );

				//				std::cout << "BACKB_DERIV1 "  <<  i << " " << atom_num_i << "     " << j << " " << atom_num_j << " " << dist_ij << " ==> " << energy << " " <<  f1(1)  << std::endl;

				F1 += f1;
				F2 += f2;
			}
		}
	} else {

		// SECOND WAY, check if this atom is a backbone oxygen atoms, cycle over other 2'-OH's.
		Size const n = find_backbone_oxygen_atom( atom_num_i );

		if ( n > 0 && rna_backbone_backbone_weight_( n ) >= 0.001 ) {

			Vector const heavy_atom_i = rsd1.xyz( atom_num_i );

			for ( graph::Graph::EdgeListConstIter
						 iter = energy_graph.get_node( i )->const_edge_list_begin();
					 iter != energy_graph.get_node( i )->const_edge_list_end();
					 ++iter ){

				Size j( ( *iter )->get_other_ind( i ) );

				if ( abs ( static_cast< int > ( i - j ) ) <= 2 ) continue;

				conformation::Residue const & rsd2( pose.residue( j ) );

				// Go over 2'-OH atoms.
				Size const m = o2prime_index_within_special_backbone_atoms_;
				Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ];

				Vector const & heavy_atom_j = rsd2.xyz( atom_num_j );
				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Distance const dist_ij = d_ij.length();

				if ( dist_ij > backbone_backbone_distance_cutoff_ ) continue;

				Real deriv( 0.0 );
				get_rna_backbone_backbone_score( dist_ij, n, deriv );

				Vector const f2( -1.0 * deriv * d_ij/ dist_ij );
				Vector const f1(  1.0 * cross( f2, heavy_atom_j ) );

				//			std::cout << "BACKB_DERIV2 "  <<  i << " " << atom_num_i << "     " << j << " " << atom_num_j << " " << dist_ij << " ==> " << energy << " " << f1(1)  << std::endl;

				F1 += f1;
				F2 += f2;


			}
		}
	}

	//	std::cout << "BACKB_DERIV_SUM " << F1( 1) << " " << F2( 1 ) << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::rna_repulsive_pair_energy(
			conformation::Residue const & rsd1,
			conformation::Residue const & rsd2
) const
{
	return ( rna_repulsive_pair_energy_one_way( rsd1, rsd2 ) +
					 rna_repulsive_pair_energy_one_way( rsd2, rsd1 ) ) ;
}

///////////////////////////////////////////////////////////////////////////////
Real
RNA_LowResolutionPotential::rna_repulsive_pair_energy_one_way(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
) const
{

	if ( !rsd1.is_RNA() ) return 0.0;
	if ( !rsd2.is_RNA() ) return 0.0;

	if ( rsd1.is_coarse() ) return 0.0;  //coarse-grained!
	if ( rsd2.is_coarse() ) return 0.0;  //coarse-grained!

	//Could use pair-moved to make this faster.

	Real rna_repulsive_score( 0.0 );

	//For speed, a cached set of atom numbers that go with RNA_backbone_oxygen_atoms_ (which is a bunch of strings)
	//rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	//	ObjexxFCL::FArray2D< Size > const &
	//		atom_numbers_for_backbone_score_calculations = rna_scoring_info.atom_numbers_for_backbone_score_calculations();

	Size const i = rsd1.seqpos();
	Size const j = rsd2.seqpos();

	//	std::string const atom_i = " OP1";
	Size const atom_num_i = atom_numbers_for_backbone_score_calculations_[ o2p_index_within_special_backbone_atoms_ ];

	Vector const heavy_atom_i = rsd1.xyz( atom_num_i );

	if ( abs( static_cast< int > ( i - j ) ) <= 2 ) return 0.0;

	// Go over sugar and phosphate oxygen atoms!
	for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ){

		if ( rna_repulsive_weight_( m ) < 0.001 ) continue;

		Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ];

		//By default only repel o2p.
		if ( !rna_repulse_all_  &&  m != o2p_index_within_special_backbone_atoms_ ) continue;

		//Don't double-count OP1 <--> OP1 interactions
		if ( atom_num_j == atom_num_i && j < i ) continue;

		Vector const heavy_atom_j( rsd2.xyz( atom_num_j ) );

		Vector const d_ij = heavy_atom_j - heavy_atom_i;

		Real const dist_ij = d_ij.length();

		if ( dist_ij > rna_repulsive_distance_cutoff_ ) continue;

		Real const score_contribution = get_rna_repulsive_score( dist_ij, m );
		rna_repulsive_score += score_contribution;

		if ( rna_verbose_  && score_contribution > 0.01 ){
			tr <<
				"REPULSIVE " <<
				rsd1.name3() <<
				I( 3, i ) << " " <<
				rsd2.atom_name( atom_num_i ) <<  " " <<
				rsd2.name3() <<
				I( 3, j ) << " " <<
				rsd2.atom_name( atom_num_j ) <<  " " <<
				F( 6, 2, score_contribution ) << " [" << dist_ij << "]" <<
				std::endl;
		}

	}

	return rna_repulsive_score;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lot of code copying here -- should minirosetta++ derivative calculation be refactored?
//   Not so much for speed but due to code redundancy in defining loops, etc.
//
void
RNA_LowResolutionPotential::eval_atom_derivative_rna_repulsive(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		Vector & F1,
 		Vector & F2 ) const
{

	Size const & i = atom_id.rsd();
	conformation::Residue const & rsd1( pose.residue( i ) );
	Size const & atom_num_i( atom_id.atomno() );

	F1 = 0.0;
	F2 = 0.0;

	if ( !rsd1.is_RNA() ) return;

	//For speed, a cached set of atom numbers that go with RNA_backbone_oxygen_atoms_ (which is a bunch of strings)
	// Need to do something different for the packer, dude. Maybe keep residue-type-specific indices
	// within the potential, determine by string lookup at the very beginning...
	//
	//	rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	//	ObjexxFCL::FArray2D< Size > const &
	//		atom_numbers_for_backbone_score_calculations = rna_scoring_info.atom_numbers_for_backbone_score_calculations();

	Size const atom_num_o2p = atom_numbers_for_backbone_score_calculations_[ o2p_index_within_special_backbone_atoms_ ];

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	if ( atom_num_i == atom_num_o2p ) {
		// FIRST WAY, check if this atom is OP1, cycle over other oxygen atoms.

		Vector const heavy_atom_i = rsd1.xyz( atom_num_i );

		for ( graph::Graph::EdgeListConstIter
					 iter = energy_graph.get_node( i )->const_edge_list_begin();
				 iter != energy_graph.get_node( i )->const_edge_list_end();
				 ++iter ){

			Size j( ( *iter )->get_other_ind( i ) );

			if ( abs( static_cast< int > ( i - j ) ) <= 2 ) continue;

			conformation::Residue const & rsd2( pose.residue( j ) );
			if ( !rsd2.is_RNA() ) continue;

			// Go over sugar and phosphate oxygen atoms!
			for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ){

				if ( rna_repulsive_weight_( m ) < 0.001 ) continue;

				//By default only repel o2p.
				if ( !rna_repulse_all_  &&  m != o2p_index_within_special_backbone_atoms_ ) continue;

				Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ];

				// Don't double-count OP1 <--> OP1 interactions
				//  if (atom_num_j == atom_num_i && j < i) continue;

				Vector const heavy_atom_j( rsd2.xyz( atom_num_j ) );

				Vector const d_ij = heavy_atom_j - heavy_atom_i;

				Real const dist_ij = d_ij.length();

				if ( dist_ij > rna_repulsive_distance_cutoff_ ) continue;

				Real deriv( 0.0 );
				get_rna_repulsive_score( dist_ij, m, deriv );

				Vector const f2( -1.0 * deriv * d_ij/ dist_ij );
				Vector const f1(  1.0 * cross( f2, heavy_atom_j ) );

				//				std::cout << "DERIV1 "  <<  i << " " << atom_num_i << "     " << j << " " << atom_num_j << " " << f1(1)  << std::endl;

				F1 += f1;
				F2 += f2;

			}

		}
	} else {

		// SECOND WAY, check if this atom is a backbone oxygen atoms, cycle over other OP1's.
		Size n = find_backbone_oxygen_atom( atom_num_i );

		if ( n > 0 && rna_repulsive_weight_( n ) >= 0.001 ) {

			Vector const heavy_atom_i = rsd1.xyz( atom_num_i );

			for ( graph::Graph::EdgeListConstIter
						 iter = energy_graph.get_node( i )->const_edge_list_begin();
					 iter != energy_graph.get_node( i )->const_edge_list_end();
					 ++iter ){

				Size j( ( *iter )->get_second_node_ind() );

				//Edges always have first node < second node. Just in case we picked the wrong one:
				if ( i == j ) j = ( *iter )->get_first_node_ind();

				if ( abs( static_cast< int > ( i - j ) ) <= 2 ) continue;

				conformation::Residue const & rsd2( pose.residue( j ) );
				if ( !rsd2.is_RNA() ) continue;

				// Go over OP1 atoms.
				Size const m = o2p_index_within_special_backbone_atoms_;
				Size const atom_num_j = atom_numbers_for_backbone_score_calculations_[ m ];

				//By default only repel o2p.
				if ( !rna_repulse_all_  &&  n != o2p_index_within_special_backbone_atoms_ ) continue;

				// Don't double-count OP1 <--> OP1 interactions
				//  if (atom_num_j == atom_num_i && j < i) continue;

				Vector const heavy_atom_j( rsd2.xyz( atom_num_j ) );

				Vector const d_ij = heavy_atom_j - heavy_atom_i;

				Real const dist_ij = d_ij.length();

				if ( dist_ij > rna_repulsive_distance_cutoff_ ) continue;

				Real deriv( 0.0 );
				get_rna_repulsive_score( dist_ij, n, deriv );

				Vector const f2( -1.0 * deriv * d_ij/ dist_ij );
				Vector const f1( cross( f2, heavy_atom_j ) );

				//				std::cout << "DERIV2 "  <<  i << " " << atom_num_i << "     " << j << " " << atom_num_j << " " << f1(1)  << std::endl;

				F1 += f1;
				F2 += f2;

			}

		}

	}

	return;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::initialize_atom_numbers_for_backbone_score_calculations()
{
	//We don't know a priori which atom numbers correspond to which
	// atom names... can we assume that all RNA residues are the same, though?

	using namespace core::chemical;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	// will Guanosine work?
	fill_atom_numbers_for_backbone_oxygens( rsd_set, na_rgu );

	//	std::cout << "HELLO?"  << check_atom_numbers_for_backbone_oxygens( rsd_set, na_rad ) << std::endl;
	//	std::cout << "HELLO? " << check_atom_numbers_for_backbone_oxygens( rsd_set, na_rcy ) << std::endl;
	//	std::cout << "HELLO? " << check_atom_numbers_for_backbone_oxygens( rsd_set, na_ura ) << std::endl;

	assert( check_atom_numbers_for_backbone_oxygens( rsd_set, na_rad ) );
	assert( check_atom_numbers_for_backbone_oxygens( rsd_set, na_rcy ) );
	assert( check_atom_numbers_for_backbone_oxygens( rsd_set, na_ura ) );

}

/////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::fill_atom_numbers_for_backbone_oxygens( chemical::ResidueTypeSetCAP & rsd_set, chemical::AA const & aa )
{
	using namespace core::chemical;

	ResidueTypeCOPs const & rsd_types( rsd_set->aa_map( aa ) );
	ResidueTypeCOP const & rsd_type = rsd_types[ 1 ]; //This better work.

	for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ) {
		atom_numbers_for_backbone_score_calculations_.push_back( rsd_type->atom_index( RNA_backbone_oxygen_atoms_[ m ] ) );
	}

}

/////////////////////////////////////////////////////
bool
RNA_LowResolutionPotential::check_atom_numbers_for_backbone_oxygens( chemical::ResidueTypeSetCAP & rsd_set, chemical::AA const & aa ) const
{
	using namespace core::chemical;

	ResidueTypeCOPs const & rsd_types( rsd_set->aa_map( aa ) );
	ResidueTypeCOP const &  rsd_type = rsd_types[ 1 ]; //This better work.

	for ( Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ) {
		if ( atom_numbers_for_backbone_score_calculations_[m] != rsd_type->atom_index( RNA_backbone_oxygen_atoms_[ m ] ) ) return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DEPRECATED.
// void
// RNA_LowResolutionPotential::initialize_atom_numbers_for_backbone_score_calculations( pose::Pose & pose ) const
// {
// 	//We don't know a priori which atom numbers correspond to which
// 	// atom names (e.g., O2' on an adenosine could be different depending
// 	// on whether its at a chainbreak, terminus, etc.)
// 	//Better to do a quick setup every time to pinpoint atoms that require
// 	//  monitoring for VDW clashes.

// 	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
// 	ObjexxFCL::FArray2D< Size > &
// 		atom_numbers_for_backbone_score_calculations( rna_scoring_info.atom_numbers_for_backbone_score_calculations() );

// 	Size const total_residue( pose.total_residue() );

// 	atom_numbers_for_backbone_score_calculations.dimension( total_residue, num_RNA_backbone_oxygen_atoms_ );
// 	atom_numbers_for_backbone_score_calculations = 0;

// 	for (Size i = 1; i <= total_residue; i++ ) {
// 		conformation::Residue const & rsd( pose.residue( i ) );

// 		if ( rsd.is_RNA() )	{
// 			for (Size m = 1; m <= num_RNA_backbone_oxygen_atoms_; m++ ) {
// 				atom_numbers_for_backbone_score_calculations( i, m ) =  rsd.atom_index( RNA_backbone_oxygen_atoms_[ m ] );
// 			}
// 		}
// 	}

// }

//////////////////////////////////////////////////////////////////////////////////
void
RNA_LowResolutionPotential::finalize( pose::Pose & ) const
{
	// ALL THIS "CALCULATED" BOOK-KEEPING WAS CONFUSING, AND NOW DOES NOT BUY ME MUCH.
// 	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );

// 	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
// 	rna_centroid_info.calculated() = false;

// 	rna::RNA_RawBaseBaseInfo & rna_raw_base_base_info( rna_scoring_info.rna_raw_base_base_info() );
// 	rna_raw_base_base_info.calculated() = false;
}

//////////////////////////////////////////////////////////
bool
RNA_LowResolutionPotential::check_clear_for_stacking(
	pose::Pose & pose,
	Size const & i,
	int const & sign
) const
{

	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//Doesn't recalculate stuff if already updated:
	rna_centroid_info.update( pose );

	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );


	//	Real const Z_CUTOFF( 2.5 );

	//	Size const total_residue = pose.total_residue();
	conformation::Residue const & res_i( pose.residue( i ) );

	if ( !res_i.is_RNA() ) return true;

	//Note that this centroid and stubs could be calculated once at the beginning of the scoring!!!
	Vector const & centroid_i( base_centroids[i] );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	for ( Size j = 1; j <= pose.total_residue(); j++ ) {
		conformation::Residue const & res_j( pose.residue( j ) );

		if ( i == j ) continue;

		for ( Size m = 1; m <= res_j.natoms(); m++ ) {

			Vector const & atom_j( res_j.xyz( m ) );

			Vector d_ij = atom_j - centroid_i;
			Real const dist_x = dot_product( d_ij, x_i );
			Real const dist_y = dot_product( d_ij, y_i );
			Real const dist_z = dot_product( d_ij, z_i );
			Real const rho2 = dist_x*dist_x + dist_y*dist_y;

			//			if (d_ij.length() < 6.0) std::cout << dist_z << " " <<  rho2 << std::endl;
			if ( ( sign * dist_z ) >= base_stack_min_height_ &&
					 ( sign * dist_z ) <= base_stack_max_height_  &&
					 rho2 < base_stack_radius2_ ) {
				//Possible BASE STACK1
				tr << "Found stacking atom on base " << i << ": " << j << " " << res_j.atom_name( m )  << std::endl;
				return false;
			} //basepair

		} //j
	} //i

	return true;
}


//////////////////////////////////////////////////////////
bool
RNA_LowResolutionPotential::check_forming_base_pair(
	pose::Pose & pose,
	Size const & i,
	Size const & j
) const
{

	rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//Doesn't recalculate stuff if already updated:
	rna_centroid_info.update( pose );

	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );


	//	Real const Z_CUTOFF( 2.5 );

	//	Size const total_residue = pose.total_residue();
	conformation::Residue const & res_i( pose.residue( i ) );
	conformation::Residue const & res_j( pose.residue( j ) );

	if ( !res_i.is_RNA() ) return false;
	if ( !res_j.is_RNA() ) return false;

	if ( i == j ) return false;

	//Note that this centroid and stubs could be calculated once at the beginning of the scoring!!!
	Vector const & centroid_i( base_centroids[i] );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	Vector const & centroid_j( base_centroids[j] );
	//kinematics::Stub const & stub_j( base_stubs[j] );

	//  Matrix const & M_j( stub_j.M );
	//	Vector const & x_j = M_j.col_x();
	//	Vector const & y_j = M_j.col_y();
	//	Vector const & z_j = M_j.col_z();
	//	Real const cos_theta = dot_product( z_i, z_j );

	Vector d_ij = centroid_j - centroid_i;
	Real const dist_x = dot_product( d_ij, x_i );
	Real const dist_y = dot_product( d_ij, y_i );
	Real const dist_z = dot_product( d_ij, z_i );
	Real const rho2 = dist_x*dist_x + dist_y*dist_y;

	//Is it a base-pair, a base-stack, or do we ignore it?
	//Size edge_bin( 1 );
	if ( ( std::abs( dist_z ) < rna_basepair_stagger_cutoff_ ) ){
		//A possible base pair
		if ( rho2 < rna_basepair_radius_cutoff2_ ){
			return true;
		}
	}


	return false;
}



} //rna
} //scoring
} //core
