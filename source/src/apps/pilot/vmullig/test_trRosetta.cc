// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_trRosetta.cc
/// @brief An integration test (with pass/fail functionality) to ensure that trRosetta is properly integrated into Rosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// protocol headers


// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>

#ifdef USE_TENSORFLOW
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/trRosetta/trRosettaProtocol_v1.hh>
#include <protocols/trRosetta/trRosettaOutputs_v1.hh>
#include <protocols/trRosetta/trRosettaOutputsBase.hh>
#include <protocols/trRosetta/trRosettaMultipleSequenceAlignment.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/file/file_sys_util.hh>
#else
#include <basic/tensorflow_manager/util.hh>
#endif //USE_TENSORFLOW

static basic::Tracer TR("test_trRosetta");


/// @brief Indicate which commandline flags are relevant to this application.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

}

#ifdef USE_TENSORFLOW

OPT_KEY( String, expected_inputs_formatted )
OPT_KEY( String, expected_distance_output_file )
OPT_KEY( String, expected_phi_output_file )
OPT_KEY( String, expected_theta_output_file )
OPT_KEY( String, expected_omega_output_file )
OPT_KEY( Integer, model_version )
OPT_KEY( Real, error_threshold)
OPT_KEY( String, input_msa_file )

/// @brief Given an input multiple sequence alignment file and its expected conversion to integers,
/// check that we're converting correctly.
bool
check_expected_inputs_v1(
	std::string const &input_msa_file,
	std::string const &expected_conversion_file
) {
	TR << "Checking conversion of " << input_msa_file << " against " << expected_conversion_file << "." << std::endl;
	// Read the expected conversions:
	int32_t max_linesize(0);
	int32_t numlines(0);
	utility::vector1< utility::vector1< int32_t > > conversions;
	utility::vector1< std::string > const conversion_lines( utility::string_split( utility::file_contents(expected_conversion_file), '\n' ) );
	for ( std::string const & line : conversion_lines ) {
		if ( line.empty() ) continue;
		std::stringstream ss( line );
		utility::vector1< int32_t > curline;
		while ( !ss.eof() ) {
			int32_t val;
			ss >> val;
			runtime_assert( !( ss.bad() || ss.fail() ) );
			curline.push_back(val);
		}
		if ( curline.empty() ) continue;
		++numlines;
		if ( curline.size() > static_cast<core::Size>(max_linesize) ) max_linesize = curline.size();

		conversions.push_back(curline);
	}

	// Parse the input msa:
	protocols::trRosetta::trRosettaMultipleSequenceAlignment msa( utility::file_contents(input_msa_file) );

	// Compare the two:
	bool success(true);
	std::string const errmsg( "Error in check_expected_inputs_v1(): ");
	utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< int32_t > > tf_vect( msa.construct_input_tensors() );
	if ( tf_vect.size() != 3 ) {
		TR.Error << errmsg << "The number of input tensors is not equal to 3!" << std::endl;
		success = false;
		return success;
	}
	if ( tf_vect[1](1) != numlines ) {
		TR.Error << errmsg << "Expected " << numlines << " sequences in the multiple sequence alignment, but got " << tf_vect[1](1) << "!" << std::endl;
		success = false;
		return success;
	}
	if ( tf_vect[2](1) != max_linesize ) {
		TR.Error << errmsg << "Expected " << max_linesize << " residues in the aligned part of the sequences, but got " << tf_vect[2](1) << "!" << std::endl;
		success = false;
		return success;
	}
	for ( core::Size iline(1); iline<=static_cast<core::Size>(numlines); ++iline ) {
		for ( core::Size jcharacter(1); jcharacter<=static_cast<core::Size>(max_linesize); ++jcharacter ) {
			if ( tf_vect[3](iline, jcharacter) != conversions[iline][jcharacter] ) {
				TR.Error << errmsg << "Expected line " << iline << " character " << jcharacter << " to be " << conversions[iline][jcharacter] << " but was actually " << tf_vect[3](iline, jcharacter) << "!" << std::endl;
				success = false;
			}
		}
	}

	if ( success ) {
		TR << "Successfully checked conversion of " << input_msa_file << " against " << expected_conversion_file << "." << std::endl;
	} else {
		TR.Error << errmsg << "Conversion of " << input_msa_file << " FAILED when checked against " << expected_conversion_file << "!" << std::endl;
	}
	return success;
}

/// @brief Parse an expected values file for model version 1.
utility::vector1< utility::vector1< core::Real > >
load_expected_vals_v1(
	std::string const & filename
) {
	utility::vector1< utility::vector1< core::Real > > outvec;

	if ( !utility::file::file_exists(filename) ) {
		return outvec; //Skip files that don't exist.
	}

	utility::vector1< std::string > const lines( utility::string_split( utility::file_contents(filename), '\n') );
	for ( std::string /*Deliberately copied.*/ line : lines ) {
		utility::strip_whitespace( line );
		if ( line.empty() || line[0] == '#' ) {
			continue;
		}
		utility::vector1< core::Real > curvec;
		utility::vector1< std::string > const linesplit( utility::split_whitespace( line ) );
		for ( std::string const & entry : linesplit ) {
			std::stringstream ss( entry );
			core::Real curval;
			ss >> curval;
			runtime_assert( !( ss.bad() || ss.fail() ) );
			curvec.push_back(curval);
		}
		outvec.push_back(curvec);
	}

	return outvec;
}

/// @brief Compare results against expected for the version 1 model.
/// @details This fuction loads data from disk!
bool
check_against_results_v1(
	core::Real const errthresh,
	std::string const &inputmsa,
	protocols::trRosetta::trRosettaOutputsBaseCOP outputs,
	std::string const &distfile,
	std::string const &phifile,
	std::string const &thetafile,
	std::string const &omegafile
) {
	using namespace protocols::trRosetta;
	std::string const errmsg( "Error in check_against_results_v1() function in test_trRosetta application: " );
	bool success(true);

	utility::fixedsizearray1< std::string const *, 4 > const expected_files{ &distfile, &phifile, &thetafile, &omegafile };
	utility::fixedsizearray1< core::Size, 4 > const expected_bins{ 37, 13, 25, 25 };
	utility::fixedsizearray1< std::string, 4 > const names { "Dist", "Phi", "Theta", "Omega" };

	//Get the input sequence length:
	core::Size min_seq_len(0);
	{
		utility::vector1< std::string > const lines( utility::string_split( utility::file_contents(inputmsa), '\n') );
		for ( std::string line /*Deliberately copied*/ : lines ) {
			utility::rstrip_whitespace(line);
			if ( line.empty() || line[0] == '>' ) { continue; }
			if ( min_seq_len==0 || line.size() < min_seq_len ) {
				min_seq_len = line.size();
			}
		}

	}
	if ( min_seq_len == 0 ) {
		success = false;
		TR.Error << errmsg + "The minimum sequence length was 0!" << std::endl;
	}

	// Check the type of the output object.
	trRosettaOutputs_v1COP v1_outputs( utility::pointer::dynamic_pointer_cast< trRosettaOutputs_v1 const >(outputs) );
	if ( v1_outputs == nullptr ) {
		success = false;
		TR.Error << errmsg + "The output object type was incorrect!" << std::endl;
		return success;
	}

	// Check each tensor in the output object.
	if ( v1_outputs->n_dist_bins() != expected_bins[1] ) {
		TR.Error << errmsg + "Found " << v1_outputs->n_dist_bins() << " dist bins instead of " << expected_bins[1] << "." << std::endl;
		success = false;
	}
	if ( v1_outputs->n_phi_bins() != expected_bins[2] ) {
		TR.Error << errmsg + "Found " << v1_outputs->n_phi_bins() << " phi bins instead of " << expected_bins[2] << "." << std::endl;
		success = false;
	}
	if ( v1_outputs->n_theta_bins() != expected_bins[3] ) {
		TR.Error << errmsg + "Found " << v1_outputs->n_theta_bins() << " theta bins instead of " << expected_bins[3] << "." << std::endl;
		success = false;
	}
	if ( v1_outputs->n_omega_bins() != expected_bins[4] ) {
		TR.Error << errmsg + "Found " << v1_outputs->n_omega_bins() << " omega bins instead of " << expected_bins[4] << "." << std::endl;
		success = false;
	}
	if ( v1_outputs->sequence_length() != min_seq_len ) {
		TR.Error << errmsg + "Found a sequence length of " << v1_outputs->sequence_length() << ", but expected " << min_seq_len << "." << std::endl;
		success = false;
	}
	bool at_least_one_file(false);
	for ( core::Size curmetric(1); curmetric <= 4; ++curmetric ) {
		for ( core::Size curbin(1); curbin <= expected_bins[curmetric]; ++curbin ) {

			//Load expected vals:
			utility::vector1< utility::vector1< core::Real > > expval( load_expected_vals_v1( *(expected_files[curmetric]) + std::to_string(curbin-1) + ".txt" ) );
			if ( expval.empty() ) continue; //Skip files that don't exist.
			at_least_one_file = true;
			if ( expval.size() != min_seq_len ) {
				TR.Error << names[curmetric] << " expected values file expected to have " << min_seq_len << " lines, but actually had " << expval.size() << "." << std::endl;
				success = false;
			}

			for ( core::Size i(1); i<=min_seq_len; ++i ) {
				if ( expval[i].size() != min_seq_len ) {
					TR.Error << "Line " << i << " of expected values file for " << names[curmetric] << " had " << expval[i].size() << " entries, but expected " << min_seq_len << "." << std::endl;
					success = false;
				}
				for ( core::Size j(1); j<=min_seq_len; ++j ) {
					core::Real curval;
					switch(curmetric) {
					case 1 :
						curval = v1_outputs->dist(i, j, curbin);
						break;
					case 2 :
						curval = v1_outputs->phi(i, j, curbin);
						break;
					case 3 :
						curval = v1_outputs->theta(i, j, curbin);
						break;
					case 4 :
						curval = v1_outputs->omega(i, j, curbin);
						break;
					default :
						utility_exit_with_message( errmsg + "Error in switch statement." );
					};
					core::Real const diff( std::abs(curval - expval[i][j] ) );
					if ( diff > errthresh ) {
						TR << names[curmetric] << " bin " << curbin << " between residues " << i << " and " << j << " was expected to be " << expval[i][j] << " but got " << curval << " instead!" << std::endl;
						success = false;
					} else {
						TR << names[curmetric] << " bin " << curbin << " between residues " << i << " and " << j << " was expected to be " << expval[i][j] << " and got " << curval << " (matches within threshold)." << std::endl;
					}
				}
			}
		}
	}

	if ( !at_least_one_file ) {
		TR.Error << errmsg + "No comparison files with the designated prefices could be found!" << std::endl;
		success = false;
	}

	return success;
}

/// @brief Read inputs, run the model, and compare outputs against expected.
bool
do_test(
	core::Size const version,
	core::Real const errthresh,
	std::string const &inputmsa,
	std::string const &distfile,
	std::string const &phifile,
	std::string const &thetafile,
	std::string const &omegafile,
	std::string const &expected_inputsfile
) {
	using namespace protocols::trRosetta;

	bool success(true);

	// Load the protocol:
	trRosettaProtocolBaseOP protocol;
	switch(version){
	case 1 :
		protocol = utility::pointer::make_shared< trRosettaProtocol_v1 >();
		break;
	default :
		TR.Error << "Error!  Bad version number " + std::to_string(version) + " provided in integration test!" << std::endl;
		success = false;
		return success;
	};

	//Read the input file:
	protocol->set_input_msa_file( inputmsa );
	trRosettaOutputsBaseCOP outputs( protocol->run() );

	// Check the results:
	switch(version){
	case 1 :
		if (
				(!check_expected_inputs_v1( inputmsa, expected_inputsfile )) ||
				(!check_against_results_v1( errthresh, inputmsa, outputs, distfile, phifile, thetafile, omegafile ) )
				) {
			success = false;
		}
		break;
	default :
		TR.Error << "Error!  Bad version number " + std::to_string(version) + " provided in integration test!" << std::endl;
		success = false;
		return success;
	};



	return success;
}

#endif //USE_TENSORFLOW

/// @brief Program entry point.
int
main(
	int argc,
	char * argv []
) {
	try {
#ifdef USE_TENSORFLOW
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::string const errmsg( "Error in test_trRosetta application: ");

		NEW_OPT( expected_inputs_formatted, "A file with the expected conversion of input sequences to integers.", "");
		NEW_OPT( expected_distance_output_file, "The prefix of the files containing expected distance output for one bin.  Filenames are of pattern prefix#.txt, where # is the 0-based bin number.", "" );
		NEW_OPT( expected_phi_output_file, "The prefix of the files containing expected phi output for one bin.  Filenames are of pattern prefix#.txt, where # is the 0-based bin number.", "" );
		NEW_OPT( expected_theta_output_file, "The prefix of the files containing expected theta output for one bin.  Filenames are of pattern prefix#.txt, where # is the 0-based bin number.", "" );
		NEW_OPT( expected_omega_output_file, "The prefix of the files containing expected omega output for one bin.  Filenames are of pattern prefix#.txt, where # is the 0-based bin number.", "" );
		NEW_OPT( model_version, "The version of the trRosetta model to use.", 1 );
		NEW_OPT( error_threshold, "The tolerance for comparing floating point numbers.  Defaults to 1.0e-4.", 1.0e-4 );
		NEW_OPT( input_msa_file, "The input multiple sequence alignment file for testing the trRosetta neural network.", "" );

		//Initialize:
		devel::init( argc, argv );

		//Get options:
		register_options();
		signed long const modelvers( option[model_version]() );
		runtime_assert_string_msg( modelvers > 0, errmsg + "The model version must be greater than zero." );
		std::string const inputmsa( option[input_msa_file]() );
		runtime_assert_string_msg( option[input_msa_file].user() && !inputmsa.empty(), errmsg + "An input multiple sequence alignment file must be provided with the -input_msa_file option." );
		std::string const expected_distfile( option[expected_distance_output_file]() );
		std::string const expected_phifile( option[expected_phi_output_file]() );
		std::string const expected_thetafile( option[expected_theta_output_file]() );
		std::string const expected_omegafile( option[expected_omega_output_file]() );
		std::string const expected_inputsfile( option[expected_inputs_formatted]() );
		runtime_assert_string_msg( option[expected_distance_output_file].user() && !expected_distfile.empty(), errmsg + "An expected distance file must be provided with the -expected_distance_output_file option." );
		runtime_assert_string_msg( option[expected_phi_output_file].user() && !expected_phifile.empty(), errmsg + "An expected phi file must be provided with the -expected_phi_output_file option." );
		runtime_assert_string_msg( option[expected_theta_output_file].user() && !expected_thetafile.empty(), errmsg + "An expected theta file must be provided with the -expected_theta_output_file option." );
		runtime_assert_string_msg( option[expected_omega_output_file].user() && !expected_omegafile.empty(), errmsg + "An expected omega file must be provided with the -expected_omega_output_file option." );
		runtime_assert_string_msg( option[expected_inputs_formatted].user() && !expected_inputsfile.empty(), errmsg + "An expected formatted inputs file must be provided with the -expected_inputs_formatted option.");
		core::Real const errthresh( option[error_threshold]() );
		runtime_assert_string_msg( errthresh > 0.0, errmsg + "The error threshold set with the -error_threshold option must be greater than zero." );

		bool const success( do_test( static_cast< core::Size >( modelvers ), errthresh, inputmsa, expected_distfile, expected_phifile, expected_thetafile, expected_omegafile, expected_inputsfile ) );
		runtime_assert_string_msg( success, errmsg + "The test_trRosetta integration test FAILED." );
		TR << "The test_trRosetta integration test PASSED." << std::endl;
#else
		//Initialize:
		devel::init( argc, argv );

		utility_exit_with_message(
			"The test_trRosetta application requires compilation with the extras=tensorflow or extras=tensorflow_gpu option.\n\n"
			+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( "test_trRosetta application" )
		);
#endif
	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
