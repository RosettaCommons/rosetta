// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   apps/pilot/ronj/pmut_scan_parallel.cc
/// @brief  Main function for running the point mutant scan protocol
/// @author Ron Jacak (ron.jacak@gmail.com)

/// MPI
#ifdef USEMPI
#include <mpi.h>
#endif

/// Core headers

// Protocol headers
#include <protocols/pmut_scan/PointMutScanDriver.hh>
#include <protocols/init.hh>

#include <basic/options/util.hh>
#include <basic/Tracer.hh>
// Auto-header: duplicate removed #include <basic/options/util.hh>

#include <utility/file/file_sys_util.hh>

#include <utility/vector1.hh>


static basic::Tracer TR("apps.pmut_scan_parallel");

// application specific options
namespace protocols {
namespace pmut_scan {

	basic::options::BooleanOptionKey const double_mutant_scan( "protocols::pmut_scan::double_mutant_scan" );
	basic::options::FileOptionKey const mutants_list( "protocols::pmut_scan::mutants_list" );
	basic::options::BooleanOptionKey const output_mutant_structures( "protocols::pmut_scan::output_mutant_structures" );

} // end namespace pmut_scan
} // end namespace protocols


int
main( int argc, char * argv [] ) {

	using namespace basic::options;

#ifdef USEMPI
	MPI_Init(&argc, &argv);
#endif

	// add application specific options to options system
	option.add( protocols::pmut_scan::double_mutant_scan, "Scan for double mutants." ).def( false );
	option.add( protocols::pmut_scan::mutants_list, "List of specific (single, double, or higher order) mutants to make." );
	option.add( protocols::pmut_scan::output_mutant_structures, "Output PDB files for the mutant poses. Default: false" ).def( false );
	
	protocols::init( argc, argv );

	//
	// concatenate -s and -l flags together to get total list of PDB files
	// The called function will die with a useful error message if neither -s or -l is specified.
	// Check to make sure all of the files exist here, too.
	//
	utility::vector1< std::string > pdb_file_names = basic::options::start_files();
	utility::vector1< std::string >::iterator input_pdb_filename, last_pdb;

	for ( input_pdb_filename = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); input_pdb_filename != last_pdb; ++input_pdb_filename ) {
		if ( !utility::file::file_exists( *input_pdb_filename ) ) {
			std::cerr << "Error. Input pdb " << *input_pdb_filename << " not found." << std::endl;
			utility_exit();
		}
	}

	// save values of options
	bool double_mutant_scan = option[ protocols::pmut_scan::double_mutant_scan ];
	std::string list_file;
	if ( option[ protocols::pmut_scan::mutants_list ].user() ) {
		list_file = option[ protocols::pmut_scan::mutants_list ].value();
	}
	bool output_mutant_structures = option[ protocols::pmut_scan::output_mutant_structures ];
	
	protocols::pmut_scan::PointMutScanDriver driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );
	driver.go();

#ifdef USEMPI
	MPI_Finalize();
#endif

}
