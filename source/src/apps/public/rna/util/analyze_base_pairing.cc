// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Analyze base pairing patterns in a silent file.


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/pose/rna/util.hh>
#include <basic/options/option.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>



// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh> // AUTO IWYU For modeler,
#include <protocols/stepwise/sampler/rna/RNA_KinematicCloser.fwd.hh> // AUTO IWYU For rna,
#include <core/scoring/rna/data/RNA_DMS_Potential.fwd.hh> // AUTO IWYU For data,

using namespace core::pose;
using namespace basic::options;


using namespace core;
using namespace protocols;
using namespace protocols::stepwise;
using namespace protocols::moves;
using namespace basic::options::OptionKeys;

///////////////////////////////////////////////////////////////////////////////
void
analyze_base_pair_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring::rna::data;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;
	using namespace core::pose::rna;
	using namespace protocols::stepwise::modeler;

	using namespace protocols::stepwise::sampler::rna;
	using namespace protocols::moves;
	using namespace core::id;
	using namespace protocols::stepwise::sampler;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD /*RNA*/ );

	FullModelInfoOP my_model;

	// input stream
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = utility::pointer::make_shared< SilentFilePoseInputStream >(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
			);
		} else {
			input = utility::pointer::make_shared< SilentFilePoseInputStream >( option[ in::file::silent ]() );
		}
	} else {
		input = utility::pointer::make_shared< PDBPoseInputStream >( option[ in::file::s ]() );
	}

	PoseOP native_pose( new Pose );
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_file( *native_pose, option[ in::file::native ].value() );
	}

	// Output: write silent file
	std::string const silent_file = option[ out::file::silent  ]();
	SilentFileOptions opts; // initialized from the command line
	SilentFileData silent_file_data( opts );
	while ( input->has_another_pose() ) {
		Pose pose;
		input->fill_pose( pose, *rsd_set );
		std::string tag = core::pose::tag_from_pose( pose );
		BinarySilentStruct s( opts, pose, tag );
		add_number_base_pairs( pose, s );
		if ( native_pose ) add_number_native_base_pairs( pose, *native_pose, s );
		// Write score only == true because this is just for analysis
		silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
	}
}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	analyze_base_pair_test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> " << std::endl;
		std::cout              << "              " << argv[0] << "  -in:file:silent <silent file> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( score::weights );
		option.add_relevant( in::file::s );
		option.add_relevant( in::file::silent );
		option.add_relevant( in::file::tags );
		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::input_res );
		option.add_relevant( full_model::cutpoint_open );
		option.add_relevant( score::weights );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}

