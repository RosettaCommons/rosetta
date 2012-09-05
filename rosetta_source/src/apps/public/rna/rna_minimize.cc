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
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/init.hh>

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

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, deriv_check )
OPT_KEY( Boolean, skip_coord_constraints )
OPT_KEY( Boolean, skip_o2star_trials )
OPT_KEY( Boolean, filter_lores_base_pairs )
OPT_KEY( Boolean, vary_geometry )


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
	rna_minimizer.deriv_check( option[ deriv_check ] );
	rna_minimizer.use_coordinate_constraints( !option[ skip_coord_constraints]() );
	rna_minimizer.skip_o2star_trials( option[ skip_o2star_trials] );
	rna_minimizer.vary_bond_geometry( option[ vary_geometry ] );


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
	using namespace basic::options;

	NEW_OPT( deriv_check, "Check analytical vs. numerical derivatives",false );
	NEW_OPT( skip_coord_constraints, "Skip first stage of minimize with coordinate constraints",false );
	NEW_OPT( skip_o2star_trials, "No O2* packing in minimizer",false );
	NEW_OPT( filter_lores_base_pairs, "Filter for models that satisfy structure parameters",false );
	NEW_OPT( vary_geometry, "Let bond lengths and angles vary from ideal",false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////
	protocols::viewer::viewer_main( my_main );

}
