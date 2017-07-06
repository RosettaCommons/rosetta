// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseRNA_Dinculeotide_Sampler_Util
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @details
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/rigid_body/util.hh>
#include <protocols/stepwise/sampler/rna/RNA_NucleosideStepWiseSampler.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/util.hh>
#include <protocols/rna/movers/RNA_LoopCloser.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/toolbox/rigid_body/util.hh>
#include <protocols/stepwise/sampler/rigid_body/EulerAngles.hh>

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


using namespace core;
using namespace core::chemical::rna;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.rigid_body.util" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace rigid_body {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::Real
get_max_centroid_to_atom_distance( utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, std::string const & atom_name ){

	using namespace core::conformation;

	runtime_assert( rsd_at_origin_list.size() >= 1 );

	Real max_distance = 0;
	for ( Size n = 1; n <= rsd_at_origin_list.size(); ++n ) {
		auto const & rsd_at_origin = rsd_at_origin_list[n];
		numeric::xyzVector< core::Real > const centroid = core::chemical::rna::get_rna_base_centroid( *rsd_at_origin, false ); //optimize by returning this by reference? Apr 10, 2010

		if ( !rsd_at_origin->has( atom_name ) ) continue;
		
		Real const distance = ( rsd_at_origin->xyz( atom_name ) - centroid ).length();

		if ( max_distance < distance ) max_distance = distance;
		TR.Debug << " sugar/base conformation num: " << n << " distance = " << distance << std::endl;
	}

	TR.Debug << "max_centroid_to_atom_distance for atom: " << atom_name << " base " << name_from_aa( ( *rsd_at_origin_list[1] ).aa() ) << ": " << max_distance << std::endl;

	return max_distance;
}


//////////////////////////////////////////////////////////////////////////////////////////
// hope to be able to deprecate this soon.
void
initialize_xyz_parameters( Distance & max_distance,
	Distance & max_distance_squared,
	int & centroid_bin_min,
	int & centroid_bin_max,
	utility::vector1< core::conformation::ResidueOP > const & moving_rsd_at_origin_list,
	Size const gap_size_to_anchor ){

	// Wait, this really should depend on identities of bases. Hmm.
	// Anyway will soon be replaced with a 'fixed' max_distance. -- rhiju, feb. 2014
	Distance C5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list, " C5'" );
	Distance O5_centroid_dist = get_max_centroid_to_atom_distance( moving_rsd_at_origin_list, " O3'" );
	Distance const Max_O3_to_C5_DIST = ( gap_size_to_anchor == 0 ) ? O3I_C5I_PLUS_ONE_MAX_DIST : O3I_C5I_PLUS_TWO_MAX_DIST;

	// If not user specified, give theoretical maximum dist between the two base's centroid, +1 is to be lenient
	if ( max_distance == 0.0 ) max_distance = Max_O3_to_C5_DIST + C5_centroid_dist + O5_centroid_dist + 1.0;
	max_distance_squared = max_distance * max_distance;

	centroid_bin_min = int(  - max_distance/STANDARD_CENTROID_BIN_SIZE );
	centroid_bin_max = int( max_distance/STANDARD_CENTROID_BIN_SIZE ) - 1;
}

////////////////////////////////////////////////////////////////
// hope to be able to deprecate this soon.
utility::vector1 < core::pose::PoseOP >
setup_pose_with_moving_residue_alternative_list(
	pose::Pose const & pose,
	Size const & moving_res,
	bool const extra_chi,
	bool const use_phenix_geo
){
	using namespace sampler::rna;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace core::chemical::rna;

	//FANG: Why this does not allow syn pyrimidine by option? Does it matter?
	ChiState chi_state =
		( pose.residue_type( moving_res ).is_purine() ) ? ANY_CHI: ANTI;

	PuckerState pucker_state = ANY_PUCKER;

	RNA_NucleosideStepWiseSampler sampler( moving_res, pucker_state, chi_state );
	sampler.set_extra_chi( extra_chi );
	sampler.set_idealize_coord( use_phenix_geo );
	sampler.set_skip_same_pucker( use_phenix_geo );
	sampler.init();

	utility::vector1< pose::PoseOP > pose_clone_list;
	for ( sampler.reset(); sampler.not_end(); ++sampler ) {
		PoseOP pose_clone = pose.clone();
		sampler.apply( *pose_clone );
		pose_clone_list.push_back( pose_clone );
		//  break; // DO NOT CHECK IN! Used to check performance hit with > 1 moving_rsd
	}
	return pose_clone_list;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
DOF_bin_value( std::map< BaseBin, int, compare_base_bin > ::const_iterator const & base_bin_it, std::string const & DOF ){

	if ( DOF == "x" ) {
		return base_bin_it->first.centroid_x;
	} else if ( DOF == "y" ) {
		return base_bin_it->first.centroid_y;
	} else if ( DOF == "z" ) {
		return base_bin_it->first.centroid_z;
	} else if ( DOF == "alpha" ) {
		return base_bin_it->first.euler_alpha;
	} else if ( DOF == "euler_z" ) {
		return base_bin_it->first.euler_z;
	} else if ( DOF == "gamma" ) {
		return base_bin_it->first.euler_gamma;
	} else {
		utility_exit_with_message( "Invalid DOF = " + DOF );
		exit( 1 ); //prevent compiler warning
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
DOF_bin_size( std::string const & DOF ){

	if ( DOF == "x" || DOF == "y" || DOF == "z" ) {
		return STANDARD_CENTROID_BIN_SIZE;
	} else if ( DOF == "alpha" || DOF == "gamma" ) {
		return STANDARD_EULER_ANGLE_BIN_SIZE;
	} else if ( DOF == "euler_z" ) {
		return STANDARD_EULER_Z_BIN_SIZE;
	} else {
		utility_exit_with_message( "Invalid DOF = " + DOF );
		exit( 1 ); //prevent compiler warning
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
analyze_base_bin_map( std::map< BaseBin, int, compare_base_bin > const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two, std::string const & foldername ){

	std::map< std::pair < int, int >, int, compare_int_pair > count_density_map;
	std::map< std::pair < int, int >, int, compare_int_pair > ::iterator count_density_it, cdend;

	std::map< BaseBin, int, compare_base_bin > ::const_iterator base_bin_it, end;

	int total_count = 0;
	int total_occupied_bin = 0;

	for ( base_bin_it = base_bin_map.begin(), end = base_bin_map.end(); base_bin_it != end; ++base_bin_it ) {

		total_occupied_bin++;
		total_count = total_count + base_bin_it->second;

		auto  const & DOF_pair = std::make_pair( DOF_bin_value( base_bin_it, DOF_one ), DOF_bin_value( base_bin_it, DOF_two ) );

		count_density_it = count_density_map.find( DOF_pair );

		if ( count_density_it == count_density_map.end() ) {
			count_density_map[DOF_pair] = 1;
		} else {
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

	for ( count_density_it = count_density_map.begin(), cdend = count_density_map.end(); count_density_it != cdend; ++count_density_it ) {
		int const & DOF_one_bin_value = count_density_it->first.first;
		int const & DOF_two_bin_value = count_density_it->first.second;

		if ( DOF_one_bin_value > DOF_one_bin_max ) DOF_one_bin_max = DOF_one_bin_value;
		if ( DOF_two_bin_value > DOF_two_bin_max ) DOF_two_bin_max = DOF_two_bin_value;
		if ( DOF_one_bin_value < DOF_one_bin_min ) DOF_one_bin_min = DOF_one_bin_value;
		if ( DOF_two_bin_value < DOF_two_bin_min ) DOF_two_bin_min = DOF_two_bin_value;

	}

	for ( int DOF_one_bin_value = ( DOF_one_bin_min - 5 ); DOF_one_bin_value < ( DOF_one_bin_max + 5 ); DOF_one_bin_value++ ) {
		for ( int DOF_two_bin_value = ( DOF_two_bin_min - 5 ); DOF_two_bin_value < ( DOF_two_bin_max + 5 ); DOF_two_bin_value++ ) {

			Real const DOF_one_value = DOF_one_bin_value * DOF_bin_size( DOF_one );
			Real const DOF_two_value = DOF_two_bin_value * DOF_bin_size( DOF_two );

			int occupied_bin_count;
			std::pair < int, int > const & DOF_pair = std::make_pair( DOF_one_bin_value, DOF_two_bin_value );
			count_density_it = count_density_map.find( DOF_pair );

			if ( count_density_it == count_density_map.end() ) {
				occupied_bin_count = 0;
			} else {
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
analyze_base_bin_map( std::map< BaseBin, int, compare_base_bin > const & base_bin_map, std::string const & foldername ){

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

	TR.Debug << "centroid_bin_size = " << STANDARD_CENTROID_BIN_SIZE << "  euler_angle_bin_size = " <<  STANDARD_EULER_ANGLE_BIN_SIZE << "  euler_z_bin_size = " << STANDARD_EULER_Z_BIN_SIZE << std::endl;
}


} //rigid_body
} //rna
} //modeler
} //stepwise
} //protocols
