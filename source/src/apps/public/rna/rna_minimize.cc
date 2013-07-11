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
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/init/init.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/rna/RNA_Minimizer.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/toolbox/AllowInsert.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key


///////////////////////////////////////////////////////////////////////////////
void
rna_fullatom_minimize_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace protocols::toolbox;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

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
		protocols::rna::ensure_phosphate_nomenclature_matches_mini( native_pose );
		native_exists = true;
	}

	// minimizer setup
	protocols::rna::RNA_Minimizer rna_minimizer;
	rna_minimizer.deriv_check( option[ OptionKeys::rna::deriv_check ]() );
	rna_minimizer.use_coordinate_constraints( !option[ OptionKeys::rna::skip_coord_constraints]() );
	rna_minimizer.skip_o2star_trials( option[ OptionKeys::rna::skip_o2star_trials] );
	rna_minimizer.vary_bond_geometry( option[ OptionKeys::rna::vary_geometry ] );

	// Silent file output setup
	std::string const silent_file = option[ out::file::silent  ]();
	SilentFileData silent_file_data;

	pose::Pose pose,start_pose;

	Size i( 0 );

	while ( input->has_another_pose() ){

		input->fill_pose( pose, *rsd_set );
		protocols::rna::ensure_phosphate_nomenclature_matches_mini( pose );
		i++;
		protocols::rna::figure_out_reasonable_rna_fold_tree( pose );
		core::pose::full_model_info::fill_full_model_info_from_command_line( pose ); // only does something if -in:file:fasta specified.

		std::cout << pose.fold_tree() << std::endl;

		if ( option[ in::file::minimize_res ].user() ){
			// don't allow anything to move, and then supply minimize_res as 'extra' minimize_res.
			AllowInsertOP allow_insert = new AllowInsert( pose );
			allow_insert->set( false );
			rna_minimizer.set_allow_insert( allow_insert );
			rna_minimizer.set_extra_minimize_res( option[ in::file::minimize_res ]() );
		}


		// graphics viewer.
		if ( i == 1 ) protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

		// do it
		pose::Pose pose_init = pose;
		rna_minimizer.apply( pose );

		// tag
		std::string tag = tag_from_pose( pose );
		Size pos = tag.find( ".pdb" ); 		// remove ".pdb"
		if ( pos != std::string::npos ) tag.replace( pos, 4, "" );
		tag += "_minimize";

		BinaryRNASilentStruct s( pose, tag );

		if ( native_exists ){
			Real const rmsd_init = all_atom_rmsd( native_pose, pose_init );
			Real const rmsd      = all_atom_rmsd( native_pose, pose );
			std::cout << "All atom rmsd: " << rmsd_init  << " to " << rmsd << std::endl;
			s.add_energy( "rms", rmsd );
			s.add_energy( "rms_init", rmsd_init );
		}

		std::cout << "Outputting " << tag << " to silent file: " << silent_file << std::endl;
		silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

		std::string const out_file =  tag + ".pdb";
		dump_pdb( pose, out_file );

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
	option.add_relevant( OptionKeys::rna::skip_coord_constraints );
	option.add_relevant( OptionKeys::rna::skip_o2star_trials );
	option.add_relevant( OptionKeys::rna::deriv_check );
	option.add_relevant( in::file::minimize_res );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////
	protocols::viewer::viewer_main( my_main );
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
}
}
