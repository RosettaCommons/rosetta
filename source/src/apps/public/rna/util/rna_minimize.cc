// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/RNA_SecStruct.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/FullModelPoseBuilder.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/rna/denovo/RNA_DeNovoPoseInitializer.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.hh>
#include <core/import_pose/options/RNA_MinimizerOptions.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/stepwise/modeler/util.hh> // for other_pose.
#include <protocols/stepwise/modeler/align/util.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

OPT_KEY( String,  params_file )
OPT_KEY( Boolean,  one_torsion_test )

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key

using namespace core::io::silent;
using namespace core::pose::rna;

/*
* /////////////////////////////////////////////////////////////////
std::map< std::string, bool >
initialize_tag_is_done( std::string const & silent_file ) {

std::map< std::string, bool > tag_is_done;
utility::vector1< std::string > tags_done;
SilentFileOptions opts;
SilentFileData silent_file_data( opts );

if ( utility::file::file_exists( silent_file ) ) {
tags_done = silent_file_data.read_tags_fast( silent_file );
for ( auto const & elem : tags_done ) {
tag_is_done[ elem ] = true;
}
}

return tag_is_done;
}
*/

///////////////////////////////////////////////////////////////////////////////
void
rna_fullatom_minimize_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::import_pose::pose_stream;
	using namespace core::import_pose;
	using namespace core::pose::full_model_info;
	using namespace core::pose::toolbox;
	using namespace protocols::rna::denovo;
	using namespace protocols::rna::denovo::movers;
	using namespace protocols::stepwise;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// is_dump_pdb
	bool is_dump_pdb( true );
	if ( option[ in::file::silent ].user() ) is_dump_pdb = false;
	if ( option[ out::pdb ].user() ) is_dump_pdb = option[ out::pdb ];

	// input stream
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			SilentFilePoseInputStreamOP input1( new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
				) );
			input1->set_order_by_energy( true );
			input = input1;
		} else {
			SilentFilePoseInputStreamOP input1( new SilentFilePoseInputStream( option[ in::file::silent ]() ) );
			input1->set_order_by_energy( true );
			input = input1;
		}
	} else {
		input = PoseInputStreamOP( new PDBPoseInputStream( option[ in::file::s ]() ) );
	}

	// native pose setup -- also align pose for the sake of StepWisePoseAligner-compatible setup
	pose::PoseOP native_pose = nullptr;
	pose::PoseOP align_pose = nullptr;

	using namespace core::import_pose::options;
	// minimizer setup
	RNA_MinimizerOptionsOP options( new RNA_MinimizerOptions );
	options->initialize_from_command_line();
	RNA_Minimizer rna_minimizer( options );

	// Silent file output setup
	std::string const & silent_file = option[ out::file::silent  ]();

	// Only if -overwrite
	if ( option[ out::overwrite ]() ) {
		remove_silent_file_if_it_exists( silent_file );
	}

	std::map< std::string, bool > tag_is_done;
	tag_is_done = initialize_tag_is_done( silent_file );

	SilentFileOptions opts; // initialized from the command line
	SilentFileData silent_file_data( opts );

	// other poses -- for scoring collections of poses connected by (virtual) loops, using full_model_info.
	utility::vector1< pose::PoseOP > other_poses;
	if ( option[ full_model::other_poses ].user() ) get_other_poses( other_poses, option[ full_model::other_poses ](), rsd_set );

	pose::PoseOP pose( new pose::Pose );
	pose::Pose start_pose;

	Size i( 0 );

	if ( option[params_file].user() ) utility_exit_with_message( " -params_file not supported in rna_minimize anymore." );

	while ( input->has_another_pose() ) {

		input->fill_pose( *pose, *rsd_set );
		i++;

		// tag -- get this FAST so we can skip if possible.
		std::string tag = tag_from_pose( *pose );
		Size pos = tag.find( ".pdb" );   // remove ".pdb"
		if ( pos != std::string::npos ) tag.replace( pos, 4, "" );
		tag += "_minimize";

		if ( tag_is_done.find( tag ) != tag_is_done.end() && tag_is_done[ tag ] ) {
			continue;
		}

		if ( !option[ in::file::silent ].user() ) {
			cleanup( *pose );

			utility::vector1< pose::PoseOP > input_poses;
			input_poses.push_back( pose );
			FullModelPoseBuilder builder;
			builder.set_input_poses( input_poses );
			builder.set_options( option );
			builder.initialize_further_from_options();
			builder.build(); // should update input_poses[1] which is pose
		}

		// We only want to do this at most once, but we must do it once
		// pose exists. Doesn't have to be the EXACT one because in this
		// use case I think they should all have the same FMI.

		// In theory, if you are trying to minimize multiple silent files
		// of different modeling problems in one executable, this is trouble.
		if ( !native_pose && option[ in::file::native ].user() ) {
			core::import_pose::initialize_native_and_align_pose( native_pose, align_pose, rsd_set, pose );
		}

		if ( option[OptionKeys::constraints::cst_fa_file].user() ) {
			// Just Reads the first cst_file...
			// Not sure why but the constraint depends on the start pose given.
			// Initialize a new pose to avoid the instability.
			core::pose::Pose test_pose;
			core::pose::make_pose_from_sequence( test_pose, pose->annotated_sequence(), *rsd_set );
			ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints(
				option[OptionKeys::constraints::cst_fa_file][1], ConstraintSetOP( new ConstraintSet ), test_pose );
			pose->constraint_set( cst_set );
		}

		AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( *pose ) );
		if ( option[ in::file::minimize_res ].user() ) {
			// don't allow anything to move, and then supply minimize_res as 'extra' minimize_res.
			atom_level_domain_map->set( false );
			options->set_extra_minimize_res( option[ in::file::minimize_res ]() );
			rna_minimizer.set_options( options );
		} else if ( option[ one_torsion_test ]() ) {
			// atom_level_domain_map->set( false );
			// for ( Size n = 1; n <= pose.residue(1).natoms(); n++ ) atom_level_domain_map->set_domain( id::AtomID( n, 1 ), 1 );
			// for ( Size n = 1; n <= pose.residue(2).natoms(); n++ ) atom_level_domain_map->set_domain( id::AtomID( n, 2 ), 2 );
			// atom_level_domain_map->set( id::NamedAtomID( " P  ", 2 ), pose, false );
			// atom_level_domain_map->set( id::NamedAtomID( " O5'", 2 ), pose, false );
		} else {
			atom_level_domain_map->set( true );
		}
		rna_minimizer.set_atom_level_domain_map( atom_level_domain_map );

		// Base pair constraints
		if ( option[ OptionKeys::rna::denovo::secstruct ].user() || option[  OptionKeys::rna::denovo::secstruct_file ].user() ) {
			// Require fasta
			FullModelParametersOP full_model_parameters;
			vector1< Size > cutpoint_open_in_full_model;
			full_model_parameters = core::import_pose::get_sequence_information( option[ OptionKeys::in::file::fasta ].value()[1], cutpoint_open_in_full_model, false );

			std::string sequence = full_model_parameters->full_sequence();
			RNA_SecStruct secstruct( option[ OptionKeys::rna::denovo::secstruct ](), option[ OptionKeys::rna::denovo::secstruct_file ](), sequence );
			auto pairings = secstruct.base_pairs();
			core::pose::rna::setup_base_pair_constraints( *pose, pairings );
		}

		// graphics viewer.
		if ( i == 1 ) protocols::viewer::add_conformation_viewer( pose->conformation(), "current", 400, 400 );

		// do it
		pose::Pose pose_init = *pose;
		rna_minimizer.apply( *pose );

		BinarySilentStruct s( opts, *pose, tag );

		if ( native_pose ) {
			Real rmsd;

			rmsd = protocols::stepwise::modeler::align::superimpose_with_stepwise_aligner( *pose, *native_pose, option[ OptionKeys::stepwise::superimpose_over_all ]() );
			std::cout << "All atom rmsd over moving residues: " << tag << " " << rmsd << std::endl;
			s.add_energy( "rms", rmsd );
		}

		std::cout << "Outputting " << tag << " to silent file: " << silent_file << std::endl;
		silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

		std::string const out_file =  tag + ".pdb";
		if ( is_dump_pdb ) {
			dump_pdb( *pose, out_file );
		}
	}
}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_fullatom_minimize_test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> " << std::endl;
		std::cout              << "              " << argv[0] << "  -in:file:silent <silent file> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		utility::vector1< Size > blank_size_vector;

		option.add_relevant( OptionKeys::rna::vary_geometry );
		option.add_relevant( OptionKeys::rna::denovo::minimize::skip_coord_constraints );
		option.add_relevant( OptionKeys::rna::denovo::minimize::skip_o2prime_trials );
		option.add_relevant( OptionKeys::rna::denovo::minimize::deriv_check );
		option.add_relevant( OptionKeys::rna::denovo::minimize::min_type );
		option.add_relevant( OptionKeys::rna::denovo::minimize::minimize_use_nblist );
		option.add_relevant( OptionKeys::rna::denovo::minimize::minimize_nblist_auto_update );
		option.add_relevant( OptionKeys::rna::denovo::minimize::minimize_max_iter );
		option.add_relevant( OptionKeys::constraints::cst_fa_file );
		option.add_relevant( in::file::minimize_res );
		option.add_relevant( score::weights );
		option.add_relevant( out::pdb );
		NEW_OPT( params_file, "Input file for pairings", "" );
		NEW_OPT( one_torsion_test, "tracking down problem with geom_sol", false );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

