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

	runtime_assert( rsd.is_RNA() );

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
get_euler_stub_bin( numeric::xyzVector< core::Real > const & centroid, rotamer_sampler::rigid_body::EulerAngles const & euler_angles ){

	Base_bin base_bin;
	base_bin.centroid_x = int( centroid[0]/centroid_bin_size );
	base_bin.centroid_y = int( centroid[1]/centroid_bin_size );
	base_bin.centroid_z = int( centroid[2]/centroid_bin_size );

	base_bin.euler_alpha = int( euler_angles.alpha()/euler_angle_bin_size );
	base_bin.euler_gamma = int( euler_angles.gamma()/euler_angle_bin_size );
	base_bin.euler_z = int( euler_angles.z()/euler_z_bin_size );

	if ( centroid[0] < 0 ) base_bin.centroid_x--;
	if ( centroid[1] < 0 ) base_bin.centroid_y--;
	if ( centroid[2] < 0 ) base_bin.centroid_z--;
	if ( euler_angles.alpha() < 0 ) base_bin.euler_alpha--;
	if ( euler_angles.gamma() < 0 ) base_bin.euler_gamma--;
	if ( euler_angles.z() < 0 ) base_bin.euler_z--;

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

	TR.Debug << "centroid_bin_size = " << centroid_bin_size << "  euler_angle_bin_size = " <<  euler_angle_bin_size << "  euler_z_bin_size = " << euler_z_bin_size << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
core::Real
get_max_centroid_to_atom_distance( utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, std::string const atom_name ){

	using namespace core::conformation;

	runtime_assert( rsd_at_origin_list.size() >= 1 );

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

}
}
}
