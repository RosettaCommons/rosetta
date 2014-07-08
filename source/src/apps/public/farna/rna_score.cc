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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/data/RNA_ChemicalMappingEnergy.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/stepwise/full_model_info/FullModelInfoSetupFromCommandLine.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/util.hh>
#include <core/init/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/RNA_DataReader.hh>
#include <core/pose/PDBInfo.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

OPT_KEY( String,  params_file )

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

OPT_KEY( StringVector, original_input )

///////////////////////////////////////////////////////////////////////////////
void
rna_score_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::sampling;
	using namespace protocols::stepwise::full_model_info;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD /*RNA*/ );

	// input stream
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = new SilentFilePoseInputStream(
																						option[ in::file::silent ](),
																						option[ in::file::tags ]()
																						);
		} else {
			input = new SilentFilePoseInputStream( option[ in::file::silent ]() );
		}
	} else {
		input = new PDBPoseInputStream( option[ in::file::s ]() );
	}

	// native pose setup
	pose::Pose native_pose;
	bool native_exists( false );
	if ( option[ in::file::native ].user() ) {
		std::string native_pdb_file  = option[ in::file::native ];
		core::import_pose::pose_from_pdb( native_pose, *rsd_set, native_pdb_file );
		cleanup( native_pose );
		native_exists = true;
	}


	// score function setup
	core::scoring::ScoreFunctionOP scorefxn;
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
	}

	// Silent file output setup
	std::string const silent_file = option[ out::file::silent  ]();
	SilentFileData silent_file_data;

	FullModelInfoOP my_model;
	utility::vector1< pose::PoseOP > other_poses;

	if ( option[ original_input ].user() ) {
		utility::vector1< std::string > const & original_files = option[ original_input ]();
		utility::vector1< pose::PoseOP > original_poses;

		for ( Size n = 1; n <= original_files.size(); n++ ) {
			original_poses.push_back( get_pdb_and_cleanup( original_files[ n ], rsd_set ) );
		}
		if ( option[ full_model::other_poses ].user() ) get_other_poses( original_poses, option[ full_model::other_poses ](), rsd_set );

		//FullModelInfo (minimal object needed for add/delete)
		fill_full_model_info_from_command_line( original_poses );
		my_model = const_full_model_info( *original_poses[ 1 ] ).clone_info();

	} else {
		// other poses -- for scoring collections of poses connected by (virtual) loops, using full_model_info.
		if ( option[ full_model::other_poses ].user() ) get_other_poses( other_poses, option[ full_model::other_poses ](), rsd_set );
	}

	scoring::rna::data::RNA_DataInfoOP rna_data_info_save;
	core::scoring::rna::data::RNA_ChemicalMappingEnergy rna_chemical_mapping_energy;
	pose::Pose pose,start_pose;

	Size i( 0 );

	if ( native_exists ) (*scorefxn)( native_pose );

	while ( input->has_another_pose() ){

		input->fill_pose( pose, *rsd_set );
		i++;

		protocols::farna::RNA_StructureParameters parameters;
                if ( option[params_file].user() ) {
                        parameters.initialize(
                                        pose, option[params_file],
                                        basic::database::full_name("sampling/rna/1jj2_RNA_jump_library.dat"),
                                        false /*ignore_secstruct*/
                        );
                        // parameters.set_suppress_bp_constraint( 1.0 );
                        parameters.setup_base_pair_constraints( pose );
                        //rna_minimizer.set_allow_insert( parameters.allow_insert() );
                }

		if ( !option[ in::file::silent ].user() ) cleanup( pose );

		if ( !option[ original_input ].user() ) {
			fill_full_model_info_from_command_line( pose, other_poses ); // only does something if -in:file:fasta specified.
		} else {
			utility::vector1< Size > resnum;
			core::pose::PDBInfoCOP pdb_info = pose.pdb_info();

			if ( pdb_info )	{
				//std::cout << std::endl << "PDB Info available for this pose..." << std::endl << std::endl;
				for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( pdb_info->number( n ) );
			} else {
				for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( n );
			}

			my_model->set_res_list( resnum );
			my_model->set_other_pose_list( other_poses );
			pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, my_model );
		}

		// graphics viewer.
		//if ( i == 1 ) protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

		// do it
		if ( ! option[ score::just_calc_rmsd]() && ! option[ OptionKeys::rna::data_file ].user() ){
			(*scorefxn)( pose );
		}

		// tag
		std::string tag = tag_from_pose( pose );
		BinarySilentStruct s( pose, tag );

		if ( native_exists ){
			//Real const rmsd      = all_atom_rmsd( native_pose, pose );
			Real const rmsd = superimpose_at_fixed_res_and_get_all_atom_rmsd( pose, native_pose );
			std::cout << "All atom rmsd over moving residues: " << tag << " " << rmsd << std::endl;
			s.add_energy( "rms", rmsd );

			// Real const rms_no_bulges = superimpose_at_fixed_res_and_get_all_atom_rmsd( pose, native_pose, true );
			// std::cout << "All atom rmsd over non-bulged moving residues: " << tag << " " << rms_no_bulges << std::endl;
			// s.add_energy( "non_bulge_rms", rms_no_bulges );

			// Stem RMSD
			if ( option[params_file].user() ) {
				std::list< Size > stem_residues( parameters.get_stem_residues( pose ) );
				if ( stem_residues.size() > 0 ) {
					Real const rmsd_stems = all_atom_rmsd( native_pose, pose, stem_residues );
					s.add_energy( "rms_stem", rmsd_stems );
					std::cout << "Stems rmsd: " << rmsd_stems << std::endl;
				}
			}
		}

		// for data_file, don't actually re-score, just compute rna_chem_map score for now.
		if ( option[ OptionKeys::rna::data_file ].user() ) { // clumsy to read this in from scratch each time.
			if ( i == 1 ){
				rna_data_info_save = protocols::farna::get_rna_data_info( pose, option[ OptionKeys::rna::data_file ] ).clone();
			} else {
				scoring::rna::nonconst_rna_scoring_info_from_pose( pose ).rna_data_info() = *rna_data_info_save;
			}
			pose.update_residue_neighbors();
			s.add_energy(  "rna_chem_map",       rna_chemical_mapping_energy.calculate_energy( pose, false /*use_low_res*/ ) );
			s.add_energy(  "rna_chem_map_lores", rna_chemical_mapping_energy.calculate_energy( pose , true /*use_low_res*/ ) );
		}


		std::cout << "Outputting " << tag << " to silent file: " << silent_file << std::endl;
		silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

		//std::string const out_file =  tag + ".pdb";
		//pose.dump_pdb( out_file );

		//		pose.dump_pdb( "test.pdb" );

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_score_test();

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

				utility::vector1< Size > blank_size_vector;
				utility::vector1< std::string > blank_string_vector;
				option.add_relevant( score::weights );
				option.add_relevant( in::file::s );
				option.add_relevant( in::file::silent );
				option.add_relevant( in::file::tags );
				option.add_relevant( in::file::native );
				option.add_relevant( in::file::fasta );
				option.add_relevant( in::file::input_res );
				option.add_relevant( full_model::cutpoint_open );
				option.add_relevant( score::just_calc_rmsd );
				option.add_relevant( rna::data_file );
				NEW_OPT( original_input, "If you want to rescore the poses using the original FullModelInfo from a SWM run, input those original PDBs here", blank_string_vector );
				NEW_OPT( params_file, "Input file for pairings", "" );

        ////////////////////////////////////////////////////////////////////////////
        // setup
        ////////////////////////////////////////////////////////////////////////////
        core::init::init(argc, argv);
				option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" );
				option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" );

        ////////////////////////////////////////////////////////////////////////////
        // end of setup
        ////////////////////////////////////////////////////////////////////////////
        protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
				return -1;
    }
}
