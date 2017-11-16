// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file symmetrize_new_ramamap.cc
/// @brief Takes a Ramachandran-style data file as input, and spits out a symmetric version (symmetrized by copying
/// the negative phi region of Ramachandran space).
/// @details Created 22 February 2017.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//General includes
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <numeric/angle.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Tracer:
static basic::Tracer TR( "apps.pilot.vmullig.symmetrize_new_ramamap" );

//Options (ugh -- global variables):
OPT_KEY (String, rama_file)
OPT_KEY (String, output_file)

/*****************
PROTOTYPES
******************/

/// @brief Set up the options for this pilot app.
void register_options();

/// @brief Read options from the options system, and check that user inputs are reasonable.
void get_options( std::string &input_filename, std::string &output_filename );

/// @brief Symmetrize the map and write output.
void symmetrize_map( std::string &input_filename, std::string &output_filename );

/*****************
FUNCTIONS
******************/

/// @brief Read options from the options system, and check that user inputs are reasonable.
void
get_options(
	std::string &input_filename,
	std::string &output_filename
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	runtime_assert_string_msg( option[rama_file].user(), "Error in symmetrize_new_ramamap application: an input file must be specified with the \"-rama_file <filename>\" commandline option." );
	input_filename = option[rama_file]();
	output_filename = option[output_file]();
}

/// @brief Symmetrize the map and write output.
void
symmetrize_map(
	std::string &input_filename,
	std::string &output_filename
) {
	utility::io::izstream instream;
	utility::io::ozstream outstream;

	instream.open( input_filename );
	runtime_assert_string_msg( instream.good(), "Error in symmetrize_new_ramamap application: could not open " + input_filename + " for read." );

	outstream.open( output_filename );

	utility::vector1< utility::vector1< core::Real > > negvals_list;
	utility::vector1< std::string > aa_identifiers;

	do {
		char line[2048];
		instream.getline( line, 2048 );
		if ( instream.eof() ) break;

		if ( line[0] == '#' || line[0] == '@' ) {
			outstream << line << std::endl;
			continue;
		}

		std::stringstream linestream( line );

		std::string aa_identifier;
		core::Real phi_bin, psi_bin;
		core::Real pval, eval;

		linestream >> aa_identifier >> phi_bin >> psi_bin >> pval >> eval;
		runtime_assert_string_msg( !linestream.fail(), "Error in symmetrize_new_ramamap application: could not parse line \"" + std::string(line) + "\"." );

		if ( phi_bin >= 0.0 ) continue; //Positive phi region should be ignored.

		char outtemp[2048];
		sprintf( outtemp, "%s\t%.0f\t%.0f\t%0.12g\t%0.12g", aa_identifier.c_str(), phi_bin, psi_bin, pval, eval );
		outstream << outtemp << std::endl;

		utility::vector1< core::Real > posvals(4);
		posvals[1] = -1.0*phi_bin - 10.0;
		posvals[2] = -1.0*psi_bin - 10.0;
		posvals[3] = pval;
		posvals[4] = eval;
		negvals_list.push_back(posvals);
		aa_identifiers.push_back( aa_identifier );

	} while( true );

	for ( core::Size i=1, imax=negvals_list.size(); i<=imax; ++i ) {
		char outtemp[2048];
		sprintf( outtemp, "%s\t%.0f\t%.0f\t%0.12g\t%0.12g", aa_identifiers[i].c_str(), negvals_list[i][1], negvals_list[i][2], negvals_list[i][3], negvals_list[i][4] );
		outstream << outtemp << std::endl;
	}

	instream.close();
	outstream.close();
}


/// @brief Set up the options for this pilot app.
void
register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( rama_file, "The input file to be symmetrized.  Required input.", "");
	NEW_OPT( output_file, "The output file.  Defaults to \"output.txt\" if not specified.", "output.txt");

}

/*****************
MAIN
******************/

int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try {
		register_options();
		devel::init(argc, argv);

		if ( TR.visible() ) {
			TR << "Starting symmetrize_new_ramamap.cc" << std::endl;
			TR << "Pilot app created 22 February 2017 by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory." << std::endl;
		}

		std::string input_filename, output_filename;
		get_options(input_filename, output_filename);
		symmetrize_map(input_filename, output_filename);

	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
}
