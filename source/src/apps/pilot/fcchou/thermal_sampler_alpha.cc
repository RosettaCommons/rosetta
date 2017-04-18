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


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/util.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/recces/sampler/rna/MC_RNA_Suite.hh>
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>
#include <protocols/recces/scratch/ThermalSamplingMover.hh>
#include <protocols/recces/scratch/thermal_sampler.hh>
#include <protocols/recces/util.hh>

#include <core/id/TorsionID.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

OPT_KEY( Boolean, recces_turner_mode )


using namespace core::pose;
using namespace basic::options;


using namespace core;
using namespace protocols;
using namespace protocols::stepwise;
using namespace protocols::moves;
using namespace basic::options::OptionKeys;
using namespace protocols::recces;
using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
void
thermal_sampler()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::modeler;
	using namespace protocols::recces::scratch;

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
			input = PoseInputStreamOP( new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
				) );
		} else {
			input = PoseInputStreamOP( new SilentFilePoseInputStream( option[ in::file::silent ]() ) );
		}
	} else {
		input = PoseInputStreamOP( new PDBPoseInputStream( option[ in::file::s ]() ) );
	}

	Pose pose;
	input->fill_pose( pose, *rsd_set );

	utility::vector1< pose::PoseOP > other_poses;
	if ( !option[ in::file::silent ].user() ) protocols::stepwise::setup::cleanup( pose );

	if ( !full_model_info_defined( pose ) || option[ in::file::fasta ].user() ) {
		protocols::stepwise::setup::fill_full_model_info_from_command_line( pose, other_poses ); // only does something if -in:file:fasta specified.
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );

	using namespace protocols::recces;

	ThermalSamplingMoverOP ts( new ThermalSamplingMover );
	ts->set_dumping_app( true );
	ts->set_recces_turner_mode( option[ recces_turner_mode ] );

	ts->apply( pose );
}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	thermal_sampler();

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

		NEW_OPT( recces_turner_mode, "run in recces_turner mode for some perverse reason?", "false" );

		utility::vector1< int > null_int_vector;
		utility::vector1< core::Real > null_real_vector;
		utility::vector1< Size > blank_size_vector;
		utility::vector1< std::string > blank_string_vector;
		option.add_relevant( score::weights );
		option.add_relevant( in::file::s );
		option.add_relevant( in::file::silent );
		option.add_relevant( in::file::tags );
		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::input_res );
		option.add_relevant( full_model::cutpoint_open );
		option.add_relevant( score::weights );

		//option.add_relevant( OptionKeys::recces::seq1 );
		//option.add_relevant( OptionKeys::recces::seq2 );
		option.add_relevant( OptionKeys::recces::n_cycle );
		option.add_relevant( OptionKeys::recces::temps );
		option.add_relevant( OptionKeys::recces::st_weights );
		option.add_relevant( OptionKeys::recces::out_prefix );
		//option.add_relevant( OptionKeys::recces::save_terms );
		//option.add_relevant( OptionKeys::recces::save_scores );
		option.add_relevant( OptionKeys::recces::dump_pdb );
		option.add_relevant( OptionKeys::recces::dump_silent );
		option.add_relevant( OptionKeys::recces::thermal_sampling::sample_residues );
		option.add_relevant( OptionKeys::recces::thermal_sampling::free_residues );
		option.add_relevant( OptionKeys::recces::thermal_sampling::angle_range_chi );
		option.add_relevant( OptionKeys::recces::thermal_sampling::angle_range_bb );
		option.add_relevant( OptionKeys::recces::thermal_sampling::chi_stdev );
		option.add_relevant( OptionKeys::recces::thermal_sampling::bb_stdev );
		option.add_relevant( OptionKeys::recces::thermal_sampling::standard_bb_stdev );
		option.add_relevant( OptionKeys::recces::out_torsions );
		option.add_relevant( OptionKeys::recces::dump_freq );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RNA_RESIDUE" );

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

