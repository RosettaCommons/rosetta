// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file format_converter.cc
/// @brief simple application for interconverting Rosetta recognized structure formats
/// @author Matthew O'Meara


// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/NullMover.hh>

// Platform Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>


using std::cout;
using std::endl;


int
main( int argc, char* argv [] ) {

	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		// define relevant options

		////////////////////////////////////
		//           GENERAL OPTIONS      //
		////////////////////////////////////

		// Additional places to look for files
		OPT( in::path::path );

		// path to rosetta_database
		OPT( in::path::database);

		OPT( in::ignore_unrecognized_res );
		OPT( in::ignore_waters );

		////////////////////////////////////
		//           INPUT OPTIONS        //
		////////////////////////////////////

		// residue type set for input: e.g. centroid or fa_standard
		OPT( in::file::residue_type_set );

		// PDB input
		OPT( in::file::s );     // a single file
		OPT( in::file::l );     // a list of files

		// Options for input silent file type
		OPT( in::file::silent );
		OPT( in::file::tags );
		OPT( in::file::silent_struct_type );
		OPT( in::file::silent_renumber );
		OPT( in::file::silent_read_through_errors );

		// Options for database input
		OPT( in::use_database );
		OPT( in::select_structures_from_database );
		OPT( inout::dbms::mode );
		OPT( inout::dbms::database_name );
		OPT( inout::dbms::pq_schema );
		OPT( inout::dbms::host );
		OPT( inout::dbms::user );
		OPT( inout::dbms::password );
		OPT( inout::dbms::port );

		////////////////////////////////////
		//           OUTPUT OPTIONS       //
		////////////////////////////////////

		// Overwrite existing output
		OPT( out::output );
		OPT( out::nooutput );
		OPT( out::overwrite );
		OPT( out::prefix );
		OPT( out::suffix );
		OPT( out::no_nstruct_label );
		OPT( out::file::renumber_pdb );

		OPT( out::file::residue_type_set );

		// Options for PDB output
		OPT( out::pdb_gz );
		OPT( out::pdb );

		// Options for output silent file types
		OPT( out::silent_gz );
		OPT( out::file::silent );
		OPT( out::file::silent_struct_type );
		OPT( out::file::silent_print_all_score_headers );

		// Options for database output
		OPT( out::use_database );
		// (for further database options see database input above)

		cout << endl;
		cout << endl;
		cout << endl;
		cout << " Rosetta Tool: format_converter - convert between Rosetta recognized structure formats" << endl;
		cout << " Usage:" << endl;
		cout << "   General:" << endl;
		cout << "      -database                                   Path to rosetta_database" << endl;
		cout << endl;
		cout << "   General Input:" << endl;
		cout << "      -in:path                                    Base path for input data" << endl;
		cout << "      -no_optH                                    Don't change positions of Hydrogen atoms! " << endl;
		cout << "      -ignore_unrecognized_res                    Alternatively, add a residue topology file." << endl;
		cout << "      -ignore_waters                              Ignore only the water molecules" << endl;
		cout << endl;
		cout << "   PDB input:" << endl;
		cout << "      -in:file:s *.pdb                            PDB input file name" << endl;
		cout << "      -in:file:l  list_of_pdbs                    File containing list of PDB input file names" << endl;
		cout << endl;
		cout << "   Silent input:" << endl;
		cout << "      -in:file:silent                             Silent input file name " << endl;
		cout << "      -in:file:tags                               Specify specific tags to be extracted, otherwise use all" << endl;
		cout << "      -in:file:fullatom                           for full atom structures " << endl;
		cout << "      -in:file:silent_struct_type <type>          Specify the input silent-file format e.g. " << endl;
		cout << "                                                    (protein_float, protein, rna, binary, score, binary_rna, pdb)" << endl;
		cout << "      -in:file:silent_renumber                    Renumber residues starting from 1" << endl;
		cout << "      -in:file:silent_read_through_errors         Try to salvage damaged silent files" << endl;
		cout << endl;
		cout << "   General Database:" << endl;
		cout << "      -inout:dbms:mode                         Specify database backend. default: 'sqlite3'" << endl;
		cout << "      -inout:dbms:database_name                If sqlite3 the filename for the database" << endl;
		cout << "      -inout:dbms:pq_schema                    For PostgreSQL, the schema namespace in the database to use" << endl;
		cout << "      -inout:dbms:host                         NOTE to use mysql or postgres as a backend:" << endl;
		cout << "      -inout:dbms:user                          compile with 'extras=mysql' or 'extras=postgres' and use" << endl;
		cout << "      -inout:dbms:password                      the '-inout:dbms:mode mysql' or" << endl;
		cout << "      -inout:dbms:port                          the '-inout:dbms:mode postgres' flag" << endl;
		cout << endl;
		cout << "   Database Input:" << endl;
		cout << "      -in:use_database                            Indicate that structures should be read from the given database" << endl;
		cout << "      -in:select_structures_from_database         An sql query to select which structures should be extracted. e.g.:" << endl;
		cout << "                                                  \"SELECT tag FROM structures WHERE tag = '7rsa';\"" << endl;
		cout << endl;
		cout << "   General Output:" << endl;
		cout << "      -out:nooutput                               Do not output any structures" << endl;
		cout << "      -out:no_nstruct_label                       Do not append _#### to tag for output structures" << endl;
		cout << "      -out:overwrite                              Overwrite structures if they already exist" << endl;
		cout << "      -out:path:all                               Path where to write output data" << endl;
		cout << "      -out:prefix                                 Prefix for output structures" << endl;
		cout << "      -out:suffix                                 Suffix for output structures" << endl;
		cout << "      -out:file:fullatom                          Force full-atom output" << endl;
		cout << endl;
		cout << "   PDB Output:" << endl;
		cout << "      -pdb                                        Output structures in PDB file format" << endl;
		cout << "      -pdb_gz                                     Output structures in zipped PDB file format" << endl;
		cout << endl;
		cout << "   Silent Output:" << endl;
		cout << "      -out:file:silent <filename>                 Output structures to specified silent file" << endl;
		cout << "      -out:file:silent_gz                         Output structures to specified silent file in zipped format" << endl;
		cout << "      -out:file:silent_struct_type <type>         Specify the input silent-file format e.g." << endl;
		cout << "                                                    (protein_float, protein, rna, binary, score, binary_rna, pdb)" << endl;
		cout << "      -out:file:silent_print_all_score_headers    Use if outputing structures for multiple sequences!!!" << endl;
		cout << "      -out:file:silent_preserve_H                 Preserve hydrogens in PDB silent-file format" << endl;
		cout << endl;
		cout << "   Database Output:" << endl;
		cout << "      -out:use_database                           Output structures to database (see General Database Options)" << endl;

		// options, random initialization
		devel::init( argc, argv );

		protocols::moves::NullMoverOP mover( new protocols::moves::NullMover );

		// start the job
		protocols::jd2::JobDistributor::get_instance()->go( mover );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
} // main
