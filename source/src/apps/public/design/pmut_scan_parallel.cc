// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/public/design/pmut_scan_parallel.cc
/// @brief  Main function for running the point mutant scan protocol
/// @author Ron Jacak (ron.jacak@gmail.com); Steven Lewis smlewi@gmail.com

/// Core headers

// Protocol headers
#include <protocols/pmut_scan/PointMutScanDriver.hh>
#include <protocols/pmut_scan/AlterSpecDisruptionDriver.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


static basic::Tracer TR( "apps.pmut_scan_parallel" );

//Note that namespacing these has no effect on their behavior - they have to be called un-namespaced!
// application specific options
namespace protocols {
namespace pmut_scan {

basic::options::BooleanOptionKey const double_mutant_scan( "protocols::pmut_scan::double_mutant_scan" );
basic::options::FileOptionKey const mutants_list( "protocols::pmut_scan::mutants_list" );
basic::options::BooleanOptionKey const output_mutant_structures( "protocols::pmut_scan::output_mutant_structures" );
basic::options::RealOptionKey const DDG_cutoff("protocols::pmut_scan::DDG_cutoff" );
basic::options::BooleanOptionKey const alter_spec_disruption_mode( "protocols::pmut_scan::alter_spec_disruption_mode" );

} // end namespace pmut_scan
} // end namespace protocols


int
main( int argc, char * argv [] ) {

	try {

		using namespace basic::options;
		// add application specific options to options system
		option.add( protocols::pmut_scan::double_mutant_scan, "Scan for double mutants." ).def( false );
		option.add( protocols::pmut_scan::mutants_list, "List of specific (single, double, or higher order) mutants to make." );
		option.add( protocols::pmut_scan::output_mutant_structures, "Output PDB files for the mutant poses. Default: false" ).def( false );
		option.add( protocols::pmut_scan::DDG_cutoff, "filter value for mutant scanning: do not bother printing mutants that do not improve the score by this much.  Negative = better score.  Does not interfere with output_mutant_structures.  Default: -1.0").def( -1.0 );
		option.add( protocols::pmut_scan::alter_spec_disruption_mode, "Use AlterSpecDisruption protocol instead.  Difference: assumes a two-chain system, and scans for mutations that weaken binding, as the first step of the alter_spec protocol.").def( false );
		devel::init( argc, argv );


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
		bool double_mutant_scan = option[ protocols::pmut_scan::double_mutant_scan ].value();
		std::string list_file;
		if ( option[ protocols::pmut_scan::mutants_list ].user() ) {
			list_file = option[ protocols::pmut_scan::mutants_list ].value();
		}
		bool output_mutant_structures = option[ protocols::pmut_scan::output_mutant_structures ].value();
		core::Real DDG_cutoff = option[ protocols::pmut_scan::DDG_cutoff ].value();

		if ( !option[ protocols::pmut_scan::alter_spec_disruption_mode ].value() ) {
			protocols::pmut_scan::PointMutScanDriver driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );
			driver.set_ddG_cutoff(DDG_cutoff);
			driver.go();
		} else { //yes, this is duplication, but it's not an OP'ed class, so we have to create the object in the if
			protocols::pmut_scan::AlterSpecDisruptionDriver driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );
			driver.set_ddG_cutoff(DDG_cutoff);
			driver.go();
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
