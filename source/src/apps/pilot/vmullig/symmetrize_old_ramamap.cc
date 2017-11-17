// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file symmetrize_old_ramamap.cc
/// @brief Takes a Ramachandran-style data file as input, and spits out a symmetric version (symmetrized by copying
/// the negative phi region of Ramachandran space).
/// @details Created 17 October 2016.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//General includes
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
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
static basic::Tracer TR( "apps.pilot.vmullig.symmetrize_old_ramamap" );

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

	runtime_assert_string_msg( option[rama_file].user(), "Error in symmetrize_old_ramamap application: an input file must be specified with the \"-rama_file <filename>\" commandline option." );
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
	runtime_assert_string_msg( instream.good(), "Error in symmetrize_old_ramamap application: could not open " + input_filename + " for read." );

	outstream.open( output_filename );

	do {
		char line[2048];
		instream.getline( line, 2048 );
		if ( instream.eof() ) break;

		std::stringstream linestream( line );
		core::Size aa_type, ss_type, phi_bin, psi_bin, counts;
		core::Real pval, eval;

		linestream >> aa_type >> ss_type >> phi_bin >> psi_bin >> counts >> pval >> eval;
		runtime_assert_string_msg( !linestream.fail(), "Error in symmetrize_old_ramamap application: could not parse line \"" + std::string(line) + "\"." );

		//Ensure that we're in the range [0, 360]:
		phi_bin = numeric::nonnegative_principal_angle_degrees(phi_bin);
		psi_bin = numeric::nonnegative_principal_angle_degrees(psi_bin);

		if ( phi_bin < 180 ) continue; //Positive phi region should be ignored.

		char outtemp[2048];
		sprintf( outtemp, "%lu\t%lu\t%lu\t%lu\t%lu\t%0.8f\t%0.8f", aa_type, ss_type, phi_bin, psi_bin, counts, pval, eval );
		outstream << outtemp << std::endl;
		core::Size const phiprime( static_cast< core::Size >( numeric::nonnegative_principal_angle_degrees( 360 - static_cast<signed long>(phi_bin) - 10 ) ) );
		core::Size const psiprime( static_cast< core::Size >( numeric::nonnegative_principal_angle_degrees( 360 - static_cast<signed long>(psi_bin) - 10 ) ) );
		sprintf( outtemp, "%lu\t%lu\t%lu\t%lu\t%lu\t%0.8f\t%0.8f", aa_type, ss_type, phiprime, psiprime, counts, pval, eval );
		outstream << outtemp << std::endl;

	} while( true );

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
			TR << "Starting symmetrize_old_ramamap.cc" << std::endl;
			TR << "Pilot app created 17 October 2016 by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory." << std::endl;
		}

		std::string input_filename, output_filename;
		get_options(input_filename, output_filename);
		symmetrize_map(input_filename, output_filename);

	} catch (utility::excn::Exception& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
}
