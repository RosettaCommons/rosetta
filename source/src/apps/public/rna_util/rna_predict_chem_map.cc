// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Rhiju, rhiju@stanford.edu

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple wrapper to output chemical mapping predictions, based on statistics & features stored in
//  RNA_DMS_Potential.cc; will be extended to other modifiers for which we are developing models.
//
// Output is based on likelihood, which is P( data | model ).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.hh>
#include <core/io/rna/RDAT.hh>
#include <core/init/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/farna/util.hh>
#include <protocols/viewer/viewers.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <numeric/xyz.functions.hh>

using namespace core;
//using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace core::chemical::rna;
using namespace core::scoring::rna::data;
using namespace ObjexxFCL::format; // AUTO USING NS

using utility::vector1;
using utility::tools::make_vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;


///////////////////////////////////////////////
// move this into a util when its working.
// TO DO: convert MATLAB code for confidence interval estimation  [pasted below, comments]
//  into C++ code.
//
// TO DO: Need to set up interpolation, using function(s) from numeric/interpolation/spline/
//


void
get_max_and_intervals( utility::vector1< Real > const & DMS_values,
											 utility::vector1< Real > const & logL_values,
											 Real & DMS_mean,
											 Real & DMS_std,
											 Real & DMS_maxL,
											 Real & DMS_interval_low,
											 Real & DMS_interval_high,
											 Real const /* confidence_level = 0.68 one-sigma */ ) {
	DMS_maxL = 0.0;
	DMS_mean = 0.0;
	DMS_std  = 0.0;
	DMS_interval_low = 0.0;
	DMS_interval_high = 0.0;

	Real max_logL( 0.0 );
	for ( Size n = 1; n <= logL_values.size(); n++ ){
		if ( n == 1 || logL_values[ n ] > max_logL ){
			max_logL = logL_values[ n ];
			DMS_maxL = DMS_values[ n ];
		}
	}
	if ( max_logL == 0.0 ) return; // sign that there's no data here.

	Real p_tot( 0.0 ), DMS_squared_mean( 0.0 );
	for ( Size n = 1; n <= logL_values.size(); n++ ){
		Real const p = exp( logL_values[ n ] );
		p_tot += p;
		DMS_mean         += p * DMS_values[ n ];
		DMS_squared_mean += p * DMS_values[ n ] * DMS_values[ n ];
	}
	DMS_mean         /= p_tot;
	DMS_squared_mean /= p_tot;

	DMS_std = std::sqrt( DMS_squared_mean - (DMS_mean * DMS_mean) );
}

// function  [min_ci, max_ci ] = get_confidence_interval( p, values, confidence_level );
// %
// % get_confidence_interval( p, smooth_DMS_values, confidence_level );
// %
// % Find points that mark a probability level above which integral of p is 95% of its
// % total integral.
// %
// if ~exist( 'values', 'var' ); values = [1:length(p)]; end;
// if ~exist( 'confidence_level', 'var' ) confidence_level = 0.95; end;
// assert( length(p) == length( values ) );

// p = p/sum(p);

// p_sort = sort( p );
// for i = 1:length( p_sort )
//   p_level = p_sort(i);
//   total_area(i) = sum( p( find( p >= p_level ) ) );
// end
// %clf;plot( p_sort, total_area ); pause;
// [~,idx] = unique( total_area );
// p_val = interp1( total_area(idx), p_sort(idx), confidence_level )

// i_vals = [];

// % 'root finder' -- where does probability distribution cross this contour level?
// for i = 1:length(p)-1
//   if ( ( p(i) <= p_val & p(i+1) > p_val ) | ...
//        ( p(i) >= p_val & p(i+1) < p_val ) )
//     i_vals = [i_vals, interp1( p( [i i+1] ), [i i+1],  p_val ) ];
//   end
// end

// [~, i_max ] = max( p );
// if isempty( find( i_vals < i_max ) ) i_vals = [i_vals, 1]; end;
// if isempty( find( i_vals > i_max ) ) i_vals = [i_vals, length(p)]; end;
// i_vals
// min_ci = interp1( [1:length(p)], values, min( i_vals ) );
// max_ci = interp1( [1:length(p)], values, max( i_vals ) );

///////////////////////////////////////////////
void
get_logL_DMS( pose::Pose & pose,
							vector1< vector1< Real > > & logL_values,
							vector1< Real >            & DMS_values,
							vector1< Size > & probed_res	) {

	logL_values.clear();
	DMS_values.clear();
	probed_res.clear();

	RNA_DMS_Potential & rna_dms_potential = core::scoring::ScoringManager::get_instance()->get_RNA_DMS_Potential();
	pose.update_residue_neighbors();
	rna_dms_potential.initialize( pose );

	// Input of DMS_values & interpolation not supported yet, but should be in the future.
	// Real DMS_value( 0.0 );
	// while( DMS_value <= 7.0 ){
	// 	DMS_values.push_back( DMS_value );
	// 	DMS_value += 0.05;
	// }
	DMS_values = rna_dms_potential.DMS_values();

	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		if ( !pose.residue_type( i ).is_RNA() ) continue;
		probed_res.push_back( i );
		logL_values.push_back( rna_dms_potential.get_logL_values( pose, i /*, DMS_values*/ ) );
	}
}


///////////////////////////////////////////////
void
predict_chem_map_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::import_pose;

	vector1 < std::string> pdb_files( option[ in::file::s ]() );
	std::string const file_path( option[ in::path::pdb ]( 1 ) );
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	Size count( 0 );
	std::string outfile;
	pose::Pose pose;
	Size total_residues( 0 );

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {

		std::string const pdb_file = pdb_files[n];

		std::string pdb_file_load = pdb_file;
		if ( file_path  != "./" ) pdb_file_load = file_path + '/' + pdb_file;
		pose_from_pdb( pose, *rsd_set, pdb_file_load );

		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << std::endl;
		std::cout << "Read in pose sequence: " << pose.annotated_sequence() << std::endl;

		vector1< vector1< Real > > logL_values;
		vector1< Real > DMS_values;
		vector1< Size > probed_res;
		get_logL_DMS( pose, logL_values, DMS_values, probed_res );

		///////////////////////////////////////////////////////
		std::string const logL_outfile = pdb_file + ".DMS.logL.txt";
		utility::io::ozstream logL_out( logL_outfile );
		logL_out << A( 10, "DMS" );
		for ( Size n = 1; n <= logL_values.size(); n++ ) {
			Size const i = probed_res[ n ];
			std::string restag = ObjexxFCL::string_of( pose.pdb_info()->chain( i ) ) + ":" + pose.residue( i ).name1() + ObjexxFCL::string_of(pose.pdb_info()->number( i ));
			logL_out << ' ' << A( 10, restag );
		}
		logL_out << std::endl;
		for ( Size k = 1; k <= DMS_values.size(); k++ ){
			logL_out << F( 10, 2, DMS_values[ k ] );
			for ( Size n = 1; n <= logL_values.size(); n++ ) logL_out << ' ' << F( 10, 2, logL_values[ n ][ k ] );
			logL_out << std::endl;
		}
		logL_out.close();

		///////////////////////////////////////////////////////
		std::string const maxL_outfile = pdb_file + ".DMS.txt";
		utility::io::ozstream maxL_out( maxL_outfile );
		vector1< Size > seqpos;
		vector1< Real > reactivity, reactivity_error;
		maxL_out << A(5, "Chain" ) << ' ' << A( 5, "Residue" ) <<
			' ' << A( 10, "DMS(mean)" ) <<
			' ' << A( 10, "DMS(std)" ) <<
			' ' << A( 10, "DMS(maxL)" ) <<
			' ' << A( 10, "DMS-1sig-low" ) << ' ' << A( 10, "DMS-1sig-hi" ) << std::endl;
		for ( Size n = 1; n <= logL_values.size(); n++ ) {
			Real DMS_mean, DMS_std, DMS_maxL, DMS_interval_low, DMS_interval_high;
			get_max_and_intervals( DMS_values, logL_values[n], DMS_mean, DMS_std, DMS_maxL, DMS_interval_low, DMS_interval_high, 0.68 );
			Size const i = probed_res[ n ];
			maxL_out << A(5,pose.pdb_info()->chain( i )) << ' ' <<
				I(5,pose.pdb_info()->number(i) )  << ' ' <<
				F( 10, 2, DMS_mean ) <<
				F( 10, 2, DMS_std ) <<
				F( 10, 2, DMS_maxL ) <<
				' ' << F( 10, 2, DMS_interval_low ) << ' ' << F(10,2,DMS_interval_high ) << std::endl;

			// for RDAT below.
			seqpos.push_back( pose.pdb_info()->number(i) );
			reactivity.push_back( DMS_mean );
			reactivity_error.push_back( DMS_std );
		}
		maxL_out.close();

		std::cout << "Outputted chemical mapping full log-likelihood profiles to: " << logL_outfile << std::endl;
		std::cout << "Outputted chemical mapping max-likelihood & 1-sigma (68%) contours to: " << maxL_outfile << std::endl;
		std::cout << "NOTE: 1-sigma confidence intervals have not yet been ported over here from MATLAB code." << std::endl;

		///////////////////////////////////////////////////////
		std::string const rdat_file = pdb_file + ".DMS.rdat";
		core::io::rna::RDAT rdat;
		rdat.fill_header_information( pose );
		rdat.set_seqpos( seqpos );
		rdat.set_data_annotations( make_vector1( make_vector1( std::make_pair("modifier","DMS" ) ) ) );
		rdat.set_reactivity( make_vector1( reactivity ) );
		rdat.set_reactivity_error( make_vector1( reactivity_error ) );
		rdat.output_rdat_to_file( rdat_file );

	}



}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	predict_chem_map_test();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}



///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {
		using namespace basic::options;

		core::init::init(argc, argv);

		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" );

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
