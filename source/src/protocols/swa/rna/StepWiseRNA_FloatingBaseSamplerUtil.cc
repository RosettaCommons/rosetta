// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Dinculeotide_Sampler_Util
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSamplerUtil.hh>
#include <protocols/rotamer_sampler/rna/RNA_NucleosideRotamer.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/id/TorsionID.hh>
//////////////////////////////////
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <set>
#include <numeric/conversions.hh>
#include <numeric/NumericTraits.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>


using namespace core;
using namespace core::chemical::rna;

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_FloatingBaseSamplerUtil" ) ;

namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_base_stack( core::kinematics::Stub const & moving_res_base,
							utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
							core::Real const base_axis_CUTOFF,
							core::Real const base_planarity_CUTOFF ){

	core::Real const small_offset = 0.000001; //0.0001;

	for ( Size i = 1; i <= other_residues_base_list.size(); i++ ){

		core::kinematics::Stub const & base_info = other_residues_base_list[i];
		numeric::xyzVector< Real > const other_z_vector = base_info.M.col_z();
		numeric::xyzVector< Real > const rebuild_z_vector = moving_res_base.M.col_z();

		numeric::xyzVector< Real > centroid_diff;
		subtract( moving_res_base.v, base_info.v, centroid_diff );
		Real const centroid_distance = centroid_diff.length();

		if ( centroid_distance > ( 6.3640 + small_offset ) ) continue;

		Real const base_z_offset_one = std::abs( dot( centroid_diff, other_z_vector ) );
		Real const base_z_offset_two = std::abs( dot( centroid_diff, rebuild_z_vector ) );

		if ( ( base_z_offset_one > ( 4.5000 + small_offset ) || base_z_offset_one < ( 2.5000 - small_offset ) ) && ( base_z_offset_two > ( 4.5000 + small_offset ) || base_z_offset_two < ( 2.5000 - small_offset ) ) ) continue;

		Real const base_axis_one = base_z_offset_one/centroid_distance;
		Real const base_axis_two = base_z_offset_two/centroid_distance;

		if ( base_axis_one < ( base_axis_CUTOFF - small_offset ) && base_axis_two < ( base_axis_CUTOFF - small_offset ) ) continue;

		Real const base_planarity = std::abs( dot( other_z_vector, rebuild_z_vector ) );

		if ( base_planarity < ( base_planarity_CUTOFF - small_offset ) ) continue;

		return true; //If reach this point means success!
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_base_pair( core::kinematics::Stub const & moving_res_base,
						utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
							core::Real const base_axis_CUTOFF,
						 core::Real const base_planarity_CUTOFF ){

	core::Real const small_offset = 0.000001; //0.0001;

	for ( Size i = 1; i <= other_residues_base_list.size(); i++ ){

		core::kinematics::Stub const & base_info = other_residues_base_list[i];
		numeric::xyzVector< Real > const other_z_vector = base_info.M.col_z();
		numeric::xyzVector< Real > const rebuild_z_vector = moving_res_base.M.col_z();

		numeric::xyzVector< Real > centroid_diff;
		subtract( moving_res_base.v, base_info.v, centroid_diff );

		Real const centroid_distance = centroid_diff.length();

		if ( centroid_distance < ( 5.0000 - small_offset ) || centroid_distance > ( 12.0000 + small_offset ) ) continue;

		Real const base_z_offset_one = std::abs( dot( centroid_diff, other_z_vector ) );
		Real const base_z_offset_two = std::abs( dot( centroid_diff, rebuild_z_vector ) );

		if ( base_z_offset_one > ( 3.0000 + small_offset ) && base_z_offset_two > ( 3.0000 + small_offset ) ) continue;

		Real const base_axis_one = base_z_offset_one/centroid_distance;
		Real const base_axis_two = base_z_offset_two/centroid_distance;

		if ( base_axis_one > ( base_axis_CUTOFF + small_offset ) && base_axis_two > ( base_axis_CUTOFF + small_offset ) ) continue; //This is a stronger condition compare to baze_z_off_set check

		Real const base_planarity = std::abs( dot( rebuild_z_vector, other_z_vector ) );

		if ( base_planarity < ( base_planarity_CUTOFF - small_offset )  ) continue;

		numeric::xyzVector< Real > const centroid_diff_parallel_one = dot( centroid_diff, other_z_vector )*other_z_vector; //messed with this on Jan 16, 2010 Parin S.
		numeric::xyzVector< Real > const centroid_diff_perpendicular_one = centroid_diff - centroid_diff_parallel_one;
		Real const rho_one = centroid_diff_perpendicular_one.length(); //length along xy plane

		numeric::xyzVector< Real > const centroid_diff_parallel_two = dot( centroid_diff, rebuild_z_vector )*rebuild_z_vector;
		numeric::xyzVector< Real > const centroid_diff_perpendicular_two = centroid_diff - centroid_diff_parallel_two;
		Real const rho_two = centroid_diff_perpendicular_two.length();

		if ( ( rho_one < ( 5.0000 - small_offset ) || rho_one > ( 10.0000 + small_offset ) ) && ( rho_two < ( 5.0000 - small_offset ) || rho_two > ( 10.0000 + small_offset ) ) ) continue;

		return true; //If reach this point means success!
	}

	//TR.Debug << std::endl;

	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_strong_base_stack( core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list ){

	Real const base_axis_CUTOFF = 0.9000;
	Real const base_planarity_CUTOFF = 0.9000;

	return is_base_stack( moving_res_base, other_residues_base_list, base_axis_CUTOFF, base_planarity_CUTOFF );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_medium_base_stack_and_medium_base_pair( core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list ){

	bool base_stack = is_base_stack( moving_res_base, other_residues_base_list, 0.7070 /*base_axis_CUTOFF*/, 0.7070 /*base_planarity_CUTOFF*/ );

	bool base_pair = is_base_pair( moving_res_base, other_residues_base_list, 0.5000 /*base_axis_CUTOFF*/, 0.7070 /*base_planarity_CUTOFF*/ );
	//value in Base_screener_class is 0.866 Sept 16 2010, Parin S.

	return ( base_stack && base_pair );

}


//CALLED ONLY IN the floating base mode!
bool
floating_base_centroid_screening( core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list, Size const num_nucleotides, StepWiseRNA_CountStruct & count_data, bool const allow_base_pair_only_screen ){

	if ( num_nucleotides > 2 ) utility_exit_with_message( "Error: num_nucleotides > 2!" );

	if ( num_nucleotides == 2 ){ //Dinucleotide.

		bool const strong_stack_base = is_strong_base_stack( moving_res_base, other_residues_base_list );

		if ( strong_stack_base ) count_data.base_stack_count++;

		bool const medium_base_stack_and_medium_base_pair = is_medium_base_stack_and_medium_base_pair( moving_res_base, other_residues_base_list );
		if ( medium_base_stack_and_medium_base_pair ) count_data.base_pairing_count++;

		bool strict_base_pair = false;
		if ( allow_base_pair_only_screen ){
			strict_base_pair = is_base_pair( moving_res_base, other_residues_base_list, 0.2588 /*base_axis_CUTOFF*/, 0.8660 /*base_planarity_CUTOFF*/ );
			if ( strict_base_pair ) count_data.strict_base_pairing_count++;
		}

		if ( strong_stack_base || medium_base_stack_and_medium_base_pair || ( allow_base_pair_only_screen && strict_base_pair ) ){
			count_data.pass_base_centroid_screen++;

			return true;
		}

		return false;

	} else{ //num_nucleotides==1, implement in Sept 16, 2010 Parin S.

		bool const regular_base_stack = is_base_stack( moving_res_base, other_residues_base_list, 0.707 /*base_axis_CUTOFF*/, 0.707 /*base_planarity_CUTOFF*/ );

		if ( regular_base_stack ) count_data.base_stack_count++;

		bool const regular_base_pair =  is_base_pair( moving_res_base, other_residues_base_list, 0.5 /*base_axis_CUTOFF*/, 0.866  /*base_planarity_CUTOFF*/ );

		if ( regular_base_pair ) count_data.base_pairing_count++;

		if ( regular_base_stack || regular_base_pair ){
			count_data.pass_base_centroid_screen++;
			return true;
		}

		return false;
	}

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

core::kinematics::Stub
get_sugar_stub( conformation::Residue const & rsd, bool const is_prepend, bool const verbose ){

		using namespace chemical;


	std::string const center_atom = ( is_prepend ) ? " C4'" : " C3'";
	std::string const x_axis_atom = ( is_prepend ) ? " H4'" : " C2'";
	std::string const y_axis_atom = ( is_prepend ) ? " C5'" : " H3'";

	if ( verbose ){
		TR.Debug << "get_sugar_stub function: ";
		output_boolean( "Is prepend = ", is_prepend, TR.Debug );
		TR.Debug << "  center_atom = " << center_atom << "  x_axis_atom = " << x_axis_atom << "  y_axis_atom = " << y_axis_atom << std::endl;
	}

	core::kinematics::Stub anchor_sugar_stub;

		assert( rsd.is_RNA() );

	Vector x, y, z;

	Vector const origin = rsd.xyz( center_atom );

	Vector const x_axis_coord = rsd.xyz( x_axis_atom );
	x = x_axis_coord - origin;
	x.normalize();

		Vector const y_axis_coord = rsd.xyz( y_axis_atom );
		y = y_axis_coord - origin; //not orthonormal yet...
	z = cross( x, y );  //Can cross here even though y is not orthogonal to x since the component of y that is parallel to x dissapear when crossed with x since cross(x, x)=0
	z.normalize(); //Choosen H4 and C5 atom specifically so that z will be roughly parallel with the helical axis.

		y = cross( z, x ); //Now y is orthonormal
		y.normalize(); //not necessary but doesn't hurt.

	anchor_sugar_stub.v = origin ;
	anchor_sugar_stub.M = Matrix::cols( x, y, z );

	return anchor_sugar_stub;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Base_bin
get_euler_stub_bin( numeric::xyzVector< core::Real > const & centriod, Euler_angles const & euler_angles ){

//		using namespace Bin_size;

	Base_bin base_bin;
	base_bin.centroid_x = int( centriod[0]/centroid_bin_size );
	base_bin.centroid_y = int( centriod[1]/centroid_bin_size );
	base_bin.centroid_z = int( centriod[2]/centroid_bin_size );

	base_bin.euler_alpha = int( euler_angles.alpha/euler_angle_bin_size );
	base_bin.euler_gamma = int( euler_angles.gamma/euler_angle_bin_size );
	base_bin.euler_z = int( euler_angles.z/euler_z_bin_size );


	if ( centriod[0] < 0 ) base_bin.centroid_x--;
	if ( centriod[1] < 0 ) base_bin.centroid_y--;
	if ( centriod[2] < 0 ) base_bin.centroid_z--;
	if ( euler_angles.alpha < 0 ) base_bin.euler_alpha--;
	if ( euler_angles.gamma < 0 ) base_bin.euler_gamma--;
	if ( euler_angles.z < 0 ) base_bin.euler_z--;

	return base_bin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
DOF_bin_value( std::map< Base_bin, int, compare_base_bin > ::const_iterator const & base_bin_it, std::string const & DOF ){

	if ( DOF == "x" ){
		return base_bin_it->first.centroid_x;
	} else if ( DOF == "y" ){
		return base_bin_it->first.centroid_y;
	} else if ( DOF == "z" ){
		return base_bin_it->first.centroid_z;
	} else if ( DOF == "alpha" ){
		return base_bin_it->first.euler_alpha;
	} else if ( DOF == "euler_z" ){
		return base_bin_it->first.euler_z;
	} else if ( DOF == "gamma" ){
		return base_bin_it->first.euler_gamma;
	} else{
		utility_exit_with_message( "Invalid DOF = " + DOF );
		exit( 1 ); //prevent compiler warning
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
DOF_bin_size( std::string const & DOF ){

	if ( DOF == "x" || DOF == "y" || DOF == "z" ){
		return centroid_bin_size;
	} else if ( DOF == "alpha" || DOF == "gamma" ){
		return euler_angle_bin_size;
	} else if ( DOF == "euler_z" ){
		return euler_z_bin_size;
	} else{
		utility_exit_with_message( "Invalid DOF = " + DOF );
		exit( 1 ); //prevent compiler warning
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
analyze_base_bin_map( std::map< Base_bin, int, compare_base_bin > const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two, std::string const foldername ){

	std::map< std::pair < int, int >, int, compare_int_pair > count_density_map;
	std::map< std::pair < int, int >, int, compare_int_pair > ::iterator count_density_it;

	std::map< Base_bin, int, compare_base_bin > ::const_iterator base_bin_it;

	int total_count = 0;
	int total_occupied_bin = 0;

	for ( base_bin_it = base_bin_map.begin(); base_bin_it != base_bin_map.end(); base_bin_it++ ){

		total_occupied_bin++;
		total_count = total_count + base_bin_it->second;

		std::pair < int, int > const & DOF_pair = std::make_pair( DOF_bin_value( base_bin_it, DOF_one ), DOF_bin_value( base_bin_it, DOF_two ) );

		count_density_it = count_density_map.find( DOF_pair );

		if ( count_density_it == count_density_map.end() ){
			count_density_map[DOF_pair] = 1;
		} else{
			count_density_it->second++;
		}
	}

	//////////////////////Output data/////////////////////////////////////////////////////////////////////////////
	std::ofstream outfile;
	std::string filename = foldername + "Bin_" + DOF_one + "_" + DOF_two + ".txt";
	outfile.open( filename.c_str() );
	Size const spacing = 14;

	outfile << std::setw( spacing ) << DOF_one;
	outfile << std::setw( spacing ) << DOF_two;
	outfile << std::setw( 30 ) << "occupied_bin_count";
	outfile << "\n";

	int DOF_one_bin_max = 0;
	int DOF_one_bin_min = 0;
	int DOF_two_bin_max = 0;
	int DOF_two_bin_min = 0;

	for ( count_density_it = count_density_map.begin(); count_density_it != count_density_map.end(); count_density_it++ ){
		int const & DOF_one_bin_value = count_density_it->first.first;
		int const & DOF_two_bin_value = count_density_it->first.second;

		if ( DOF_one_bin_value > DOF_one_bin_max ) DOF_one_bin_max = DOF_one_bin_value;
		if ( DOF_two_bin_value > DOF_two_bin_max ) DOF_two_bin_max = DOF_two_bin_value;
		if ( DOF_one_bin_value < DOF_one_bin_min ) DOF_one_bin_min = DOF_one_bin_value;
		if ( DOF_two_bin_value < DOF_two_bin_min ) DOF_two_bin_min = DOF_two_bin_value;

	}

	for ( int DOF_one_bin_value = ( DOF_one_bin_min - 5 ); DOF_one_bin_value < ( DOF_one_bin_max + 5 ); DOF_one_bin_value++ ){
		for ( int DOF_two_bin_value = ( DOF_two_bin_min - 5 ); DOF_two_bin_value < ( DOF_two_bin_max + 5 ); DOF_two_bin_value++ ){

			Real const DOF_one_value = DOF_one_bin_value*DOF_bin_size( DOF_one );
			Real const DOF_two_value = DOF_two_bin_value*DOF_bin_size( DOF_two );

			int occupied_bin_count;
			std::pair < int, int > const & DOF_pair = std::make_pair( DOF_one_bin_value, DOF_two_bin_value );
			count_density_it = count_density_map.find( DOF_pair );

			if ( count_density_it == count_density_map.end() ){
				occupied_bin_count = 0;
			} else{
				occupied_bin_count = count_density_it->second;
			}

			outfile << std::setw( spacing ) << DOF_one_value;
			outfile << std::setw( spacing ) << DOF_two_value;
			outfile << std::setw( spacing ) << occupied_bin_count;
			outfile << "\n";
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	TR.Debug << std::setw( 50 ) << std::left << "Analysis " + DOF_one + "_" + DOF_two;
	TR.Debug << " tot_count = " << std::setw( 15 ) << std::left << total_count << " tot_occ = " << std::setw( 15 ) << std::left << total_occupied_bin;
	TR.Debug << " tot_count/tot_occ_bin = " << std::setw( 5 ) << std::left << ( double( total_count )/double( total_occupied_bin ) ) << std::endl;


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
analyze_base_bin_map( std::map< Base_bin, int, compare_base_bin > const & base_bin_map, std::string const foldername ){

	if ( system( std::string( "rm -r " + foldername ).c_str() ) == -1 ) {
		TR.Error << "Shell process failed!" << std::endl;
	}
	if ( system( std::string( "mkdir -p " + foldername ).c_str() ) == -1 ) {
		TR.Error << "Shell process failed!" << std::endl;
	}

	// const Real DEGS_PER_RAD = 180. / numeric::NumericTraits<Real>::pi(); // Unused variable causes warning.

	analyze_base_bin_map( base_bin_map, "x", "y", foldername );
	analyze_base_bin_map( base_bin_map, "x", "z", foldername );
	analyze_base_bin_map( base_bin_map, "x", "alpha", foldername );
	analyze_base_bin_map( base_bin_map, "x", "euler_z", foldername );
	analyze_base_bin_map( base_bin_map, "x", "gamma", foldername );

	analyze_base_bin_map( base_bin_map, "y", "z", foldername );
	analyze_base_bin_map( base_bin_map, "y", "alpha", foldername );
	analyze_base_bin_map( base_bin_map, "y", "euler_z", foldername );
	analyze_base_bin_map( base_bin_map, "y", "gamma", foldername );

	analyze_base_bin_map( base_bin_map, "z", "alpha", foldername );
	analyze_base_bin_map( base_bin_map, "z", "euler_z", foldername );
	analyze_base_bin_map( base_bin_map, "z", "gamma", foldername );

	analyze_base_bin_map( base_bin_map, "alpha", "euler_z", foldername );
	analyze_base_bin_map( base_bin_map, "alpha", "gamma", foldername );

	analyze_base_bin_map( base_bin_map, "euler_z", "gamma", foldername );

//		TR << "is_dinucleotide= " << is_dinucleotide << std::endl;
	TR.Debug << "centroid_bin_size = " << centroid_bin_size << "  euler_angle_bin_size = " <<  euler_angle_bin_size << "  euler_z_bin_size = " << euler_z_bin_size << std::endl;


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Actually deleted this function since it was no longer in use......On May 2, copied back a version of this function from a Mar 28 Desktop mini src.
Euler_angles
get_euler_angles( numeric::xyzMatrix< core::Real > const & coordinate_matrix ){

	numeric::xyzMatrix< core::Real > const & M = coordinate_matrix;
							const Real DEGS_PER_RAD = 180. / numeric::NumericTraits < Real > ::pi();

	Euler_angles euler_angles;
	euler_angles.alpha = atan2( M.xz(),  - M.yz() ) * DEGS_PER_RAD;
	euler_angles.z = M.zz();
	euler_angles.gamma = atan2( M.zx(), M.zy() ) * DEGS_PER_RAD ;  //tan2(y,x)=gamma

//    float atan2 (       float y,       float x );
//principal arc tangent of y/x, in the interval [-pi,+pi] radians.

	//Implement a check here to make sure that the euler_angles are be used to recalculate the coordiante May 2, 2010....

/* Defination of rotation matrix at mathworld is actually the inverse of the rotation matrix....also did not correctly input the coefficient of arctan2
	euler_angles.alpha = (  - 1 )* atan2( M.zy(), M.zx() ) * DEGS_PER_RAD; //Is this correct?? Ambuiguity with the minus sign...return in the [-Pi, Pi] range. atan2(y-value, x-value)
	euler_angles.beta  = acos( M.zz() ) * DEGS_PER_RAD;  // rerun in the [0,Pi] range
	euler_angles.z = M.zz(); //Use z instead of beta to make space uniform
	euler_angles.gamma = atan2( M.yz(), M.xz() ) * DEGS_PER_RAD; //return in the [-Pi, Pi] range. atan2(y-value, x-value)
*/

	return euler_angles;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
convert_euler_to_coordinate_matrix( Euler_angles const & E, numeric::xyzMatrix< core::Real > & coordinate_matrix ){

	//Probably could save time if determine x and y and take cross product to determine z
	//Should determine using both ways and check for consistency

	coordinate_matrix.xx( cos( E.alpha )*cos( E.gamma ) - sin( E.alpha )*cos( E.beta )*sin( E.gamma ) );
	coordinate_matrix.xy(  - cos( E.alpha )*sin( E.gamma ) - sin( E.alpha )*cos( E.beta )*cos( E.gamma ) );
	coordinate_matrix.xz( sin( E.alpha )*sin( E.beta ) );
	coordinate_matrix.yx( sin( E.alpha )*cos( E.gamma ) + cos( E.alpha )*cos( E.beta )*sin( E.gamma ) );
	coordinate_matrix.yy(  - sin( E.alpha )*sin( E.gamma ) + cos( E.alpha )*cos( E.beta )*cos( E.gamma ) ); //Found bug on Feb 13, 2010...previously had cos(E.gamma) instead of sin(E.gamma)
	coordinate_matrix.yz(  - cos( E.alpha )*sin( E.beta ) );
	coordinate_matrix.zx( sin( E.beta ) *sin( E.gamma ) );
	coordinate_matrix.zy( sin( E.beta ) *cos( E.gamma ) );
	coordinate_matrix.zz( cos( E.beta ) );

	Real determinant = coordinate_matrix.det();
	//if(determinant>1.00000000000001 || determinant<0.99999999999999){ //Feb 12, 2012 This might lead to server-test error at R47200
	if ( determinant > 1.000001 || determinant < 0.999999 ){ //Feb 12, 2012 This might lead to server-test error at R47200
		utility_exit_with_message( "determinant != 1.00 !!!" );
	}
}

void
translate_then_rotate_pose( core::pose::Pose & pose, numeric::xyzVector< core::Real > const & vector, numeric::xyzMatrix< core::Real > const matrix, bool const verbose ){

	using namespace core::id;
	using namespace core::pose;


	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

		conformation::Residue const & rsd( pose.residue( seq_num ) );

		for ( Size at = 1; at <= rsd.natoms(); at++ ){
			std::string const & atom_name = rsd.type().atom_name( at );
			if ( verbose )	TR << "seq_num = " << seq_num << " atom_name = " << atom_name << std::endl;

			id::AtomID const id( at, seq_num );

			pose.set_xyz( id, pose.xyz( id ) + vector );
			pose.set_xyz( id, matrix * pose.xyz( id ) );
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set sets the base to the origin....
void
set_to_origin( pose::Pose & pose, Size const seq_num, bool verbose ){

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::id;

	conformation::Residue const & rsd( pose.residue( seq_num ) );

	numeric::xyzVector< core::Real > centroid = core::chemical::rna::get_rna_base_centroid( rsd, verbose );

	numeric::xyzMatrix< core::Real > base_coordinate_matrix = core::chemical::rna::get_rna_base_coordinate_system( rsd, centroid );

	numeric::xyzMatrix< core::Real > invert_coordinate_matrix = inverse( base_coordinate_matrix );

	for ( Size at = 1; at <= rsd.natoms(); at++ ){
		id::AtomID const id( at, seq_num );

		pose.set_xyz( id, pose.xyz( id ) - centroid ); //I think the order here does matter. Translate centroid to origin.
		pose.set_xyz( id, invert_coordinate_matrix * pose.xyz( id ) ); //I think the order here does matter. Rotate coordinate so that it equal to Roseeta internal reference frame

	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Find the " O3'" atom coordinate (prepend), " C5'" atom coordinate (append)
void
get_specific_atom_coordinate( std::string const & atom_name,
														 numeric::xyzVector< core::Real > & atom_pos,
														 core::conformation::Residue const & rsd_at_origin,
														 core::kinematics::Stub const & moving_res_base_stub ){

	numeric::xyzVector< core::Real > const & new_centroid = moving_res_base_stub.v;
	numeric::xyzMatrix< core::Real > const & new_coordinate_matrix = moving_res_base_stub.M;

	//STILL NEED TO CHECK THAT THIS NEW IMPLEMENTATION WORKS!! Apr 10, 2010 Parin
	atom_pos = new_coordinate_matrix * rsd_at_origin.xyz( atom_name ); //I think the order here does matter.
	atom_pos = atom_pos + new_centroid; //I think the order here does matter.

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::Real
get_max_centroid_to_atom_distance( utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, std::string const atom_name ){

	using namespace core::conformation;

	if ( rsd_at_origin_list.size() < 1 ){
		utility_exit_with_message( "rsd_at_origin_list.size() < 1!!" );
	}

	Real max_distance = 0;

	for ( Size n = 1; n <= rsd_at_origin_list.size(); n++ ){
		Residue const & rsd_at_origin = ( *rsd_at_origin_list[n] );
		numeric::xyzVector< core::Real > const centroid = core::chemical::rna::get_rna_base_centroid( rsd_at_origin, false ); //optimize by returning this by reference? Apr 10, 2010

		Real const distance = ( rsd_at_origin.xyz( atom_name ) - centroid ).length();

		if ( max_distance < distance ) max_distance = distance;
		TR.Debug << " sugar/base conformation num: " << n << " distance = " << distance << std::endl;
	}

	TR.Debug << "max_centroid_to_atom_distance for atom: " << atom_name << " base " << name_from_aa( ( *rsd_at_origin_list[1] ).aa() ) << ": " << max_distance << std::endl;

	return max_distance;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1 < core::conformation::ResidueOP >
setup_residue_at_origin_list(
	pose::Pose const & pose,
	Size const & moving_res,
	bool const extra_chi,
	bool const use_phenix_geo
){
	using namespace rotamer_sampler::rna;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace core::chemical::rna;

	//FANG: Why this does not allow syn pyrimidine by option? Does it matter?
	core::Size base_state =
			( is_purine( pose.residue( moving_res ) ) ) ? WHATEVER: ANTI;

	core::Size pucker_state = WHATEVER;

	RNA_NucleosideRotamer sampler( moving_res, pucker_state, base_state );
	sampler.set_extra_chi( extra_chi );
	sampler.set_idealize_coord( use_phenix_geo );
	sampler.set_skip_same_pucker( use_phenix_geo );
	sampler.init();

	utility::vector1< ResidueOP > rsd_at_origin_list;

	Pose pose_at_origin( pose );

	for ( sampler.reset(); sampler.not_end(); ++sampler ) {
		sampler.apply( pose_at_origin );
		set_to_origin( pose_at_origin, moving_res, false );
		ResidueOP rsd_at_origin = pose_at_origin.residue( moving_res ).clone();
		rsd_at_origin_list.push_back( rsd_at_origin );
	}

	return rsd_at_origin_list;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
check_floating_base_chain_closable( core::Size const & reference_res,
																		utility::vector1< core::pose::PoseOP > pose_data_list,
																		utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
																		core::kinematics::Stub const & moving_res_base_stub,
																		bool const is_prepend,
																		core::Size const gap_size ){

	using namespace core::conformation;

	for ( Size n = 1; n <= rsd_at_origin_list.size(); n++ ){

		Residue const & rsd_at_origin = ( *rsd_at_origin_list[n] );

		std::string const moving_atom_name = ( is_prepend ) ? "O3'" : " C5'";
		std::string const reference_atom_name = ( is_prepend ) ? " C5'" : "O3'";

		numeric::xyzVector< core::Real > atom_coordinate;

		get_specific_atom_coordinate( moving_atom_name, atom_coordinate, rsd_at_origin, moving_res_base_stub );

		for ( Size sugar_ID = 1; sugar_ID <= pose_data_list.size(); sugar_ID++ ){

			pose::Pose const & pose = ( *pose_data_list[sugar_ID] );
			if ( check_chain_closable( atom_coordinate, pose.residue( reference_res ).xyz( reference_atom_name ), gap_size ) ) return true;

		}

	}

	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Could make this call the pose_data_list version (right above)
bool
check_floating_base_chain_closable( core::Size const & reference_res,
																	 core::pose::Pose const & pose,
																	 utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, //this one correspond to the moving_base
																	 core::kinematics::Stub const & moving_res_base_stub,
																	 bool const is_prepend,
																	 Size const gap_size ){

	using namespace core::conformation;

	for ( Size n = 1; n <= rsd_at_origin_list.size(); n++ ){

		Residue const & rsd_at_origin = ( *rsd_at_origin_list[n] );

		std::string const moving_atom_name = ( is_prepend ) ? "O3'" : " C5'";
		std::string const reference_atom_name = ( is_prepend ) ? " C5'" : "O3'";

		numeric::xyzVector< core::Real > atom_coordinate;

		get_specific_atom_coordinate( moving_atom_name, atom_coordinate, rsd_at_origin, moving_res_base_stub );

		if ( check_chain_closable( atom_coordinate, pose.residue( reference_res ).xyz( reference_atom_name ), gap_size ) ) return true;
	}

	return false;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
set_base_coordinate_frame( pose::Pose & pose, Size const & seq_num, core::conformation::Residue const & rsd_at_origin, core::kinematics::Stub const & moving_res_base_stub ){

	utility::vector1< std::pair < id::AtomID, numeric::xyzVector< core::Real > > > xyz_list;
	get_atom_coordinates( xyz_list, seq_num, rsd_at_origin, moving_res_base_stub );

	for ( Size n = 1; n <= xyz_list.size(); n++ ){
		pose.set_xyz( xyz_list[n].first, xyz_list[n].second );
	}

}


}
}
}
