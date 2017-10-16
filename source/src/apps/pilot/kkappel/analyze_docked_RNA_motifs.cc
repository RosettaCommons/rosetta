// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Analyze RNA motifs docked onto proteins to try to design protein binding structures


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
#include <core/pose/extra_pose_info_util.hh>
#include <core/init/init.hh>
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
#include <protocols/rna/denovo/setup/RNA_DeNovoPoseInitializer.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/EtableEnergyCreator.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairGeneric.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <basic/Tracer.hh>

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
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

OPT_KEY( Boolean, native_RNA_binding_res )
OPT_KEY( IntegerVector, RNA_binding_res )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kkappel.analyze_docked_RNA_motifs" );

///////////////////////////////////////////////////////////////////////////////
void get_contacts( core::pose::Pose const & pose, std::string const & tag ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;

	// Loop through the RNA binding residues
	// count how many of them are in contact with RNA

	// contact defined as: at least one RNA atom within 3.5A of one of the protein atoms

	utility::vector1< core::Size > RNA_binding_res_in_contact;
	bool protein_res_contacted;

	// I'm sure there's a smarter way to do this
	// for each RNA binding residue
	for ( core::Size i = 1; i <= option[ RNA_binding_res ]().size(); ++i ) {
		protein_res_contacted = false;
		core::Size prot_res = option[ RNA_binding_res ]()[i];
		// loop through all the protein atoms in this residue
		for ( core::Size prot_atom = 1; prot_atom <= pose.residue( prot_res ).natoms(); ++prot_atom ) {
			Vector protein_xyz( pose.residue( prot_res ).xyz( prot_atom ));
			// loop through all the RNA residues
			for ( core::Size RNA_res = 1; RNA_res <= pose.total_residue(); ++RNA_res ) {
				if ( !pose.residue( RNA_res ).is_RNA() ) continue;
				// loop through all the atoms in the RNA residue
				for ( core::Size RNA_atom = 1; RNA_atom <= pose.residue( RNA_res ).natoms(); ++RNA_atom ) {

					Vector RNA_xyz( pose.residue( RNA_res ).xyz( RNA_atom ));
					// if the distance is greater than 20A, get out of this for loop, nothing is ever
					// going to contact the protein
					Vector r_vector = protein_xyz - RNA_xyz;
					core::Real distance( r_vector.length() );
					if ( distance > 20.0 ) break;
					if ( distance <= 3.5 ) {
						RNA_binding_res_in_contact.push_back( i );
						protein_res_contacted = true;
						break; // only need 1 for it to count as in contact
					}
				}
				if ( protein_res_contacted ) break;
			}
			if ( protein_res_contacted ) break;
		}
	}

	core::Real num_contacts = RNA_binding_res_in_contact.size();
	TR << TR.Blue << tag << " contacts " << num_contacts << "/" << option[ RNA_binding_res ]().size() << " predicted RNA binding residues" << std::endl;

}
///////////////////////////////////////////////////////////////////////////////
void get_RNA_binding_residues( core::pose::Pose const & pose, std::string const & tag ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;

	// Loop through the RNA binding residues
	// count how many of them are in contact with RNA

	// contact defined as: at least one RNA atom within 3.5A of one of the protein atoms

	utility::vector1< core::Size > protein_RNA_binding_res;
	bool protein_res_contacted;

	for ( core::Size prot_res = 1; prot_res <= pose.total_residue(); ++prot_res ) {
		if ( !pose.residue( prot_res ).is_protein() ) continue;
		protein_res_contacted = false;
		// loop through all the protein atoms in this residue
		for ( core::Size prot_atom = 1; prot_atom <= pose.residue( prot_res ).natoms(); ++prot_atom ) {
			Vector protein_xyz( pose.residue( prot_res ).xyz( prot_atom ));
			// loop through all the RNA residues
			for ( core::Size RNA_res = 1; RNA_res <= pose.total_residue(); ++RNA_res ) {
				if ( !pose.residue( RNA_res ).is_RNA() ) continue;
				// loop through all the atoms in the RNA residue
				for ( core::Size RNA_atom = 1; RNA_atom <= pose.residue( RNA_res ).natoms(); ++RNA_atom ) {

					Vector RNA_xyz( pose.residue( RNA_res ).xyz( RNA_atom ));
					// if the distance is greater than 20A, get out of this for loop, nothing is ever
					// going to contact the protein
					Vector r_vector = protein_xyz - RNA_xyz;
					core::Real distance( r_vector.length() );
					if ( distance > 20.0 ) break;
					if ( distance <= 3.5 ) {
						protein_RNA_binding_res.push_back( prot_res );
						protein_res_contacted = true;
						break; // only need 1 for it to count as in contact
					}
				}
				if ( protein_res_contacted ) break;
			}
			if ( protein_res_contacted ) break;
		}
	}

	TR << TR.Blue << tag << " RNA contacts protein residues: " << protein_RNA_binding_res << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
void
count_contacts()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD /*RNA*/ );

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

	// score function setup
	core::scoring::ScoreFunctionOP scorefxn;
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
	}

	pose::Pose pose,start_pose;

	while ( input->has_another_pose() ) {

		input->fill_pose( pose, *rsd_set );

		// tag
		std::string tag = tag_from_pose( pose );

		if ( option[ native_RNA_binding_res ]() ) {
			get_RNA_binding_residues( pose, tag);
		} else {
			get_contacts( pose, tag );
		}

		//  // Loop through the RNA binding residues
		//  // count how many of them are in contact with RNA
		//
		//  // contact defined as: at least one RNA atom within 3.5A of one of the protein atoms
		//
		//  utility::vector1< core::Size > RNA_binding_res_in_contact;
		//  bool protein_res_contacted;
		//
		//  // I'm sure there's a smarter way to do this
		//  // for each RNA binding residue
		//  for ( core::Size i = 1; i <= option[ RNA_binding_res ]().size(); ++i ) {
		//   protein_res_contacted = false;
		//   core::Size prot_res = option[ RNA_binding_res ]()[i];
		//   // loop through all the protein atoms in this residue
		//   for ( core::Size prot_atom = 1; prot_atom <= pose.residue( prot_res ).natoms(); ++prot_atom ) {
		//    Vector protein_xyz( pose.residue( prot_res ).xyz( prot_atom ));
		//    // loop through all the RNA residues
		//    for ( core::Size RNA_res = 1; RNA_res <= pose.total_residue(); ++RNA_res ) {
		//     if ( !pose.residue( RNA_res ).is_RNA() ) continue;
		//     // loop through all the atoms in the RNA residue
		//     for ( core::Size RNA_atom = 1; RNA_atom <= pose.residue( RNA_res ).natoms(); ++RNA_atom ) {
		//
		//      Vector RNA_xyz( pose.residue( RNA_res ).xyz( RNA_atom ));
		//      // if the distance is greater than 20A, get out of this for loop, nothing is ever
		//      // going to contact the protein
		//      Vector r_vector = protein_xyz - RNA_xyz;
		//      core::Real distance( r_vector.length() );
		//      if ( distance > 20.0 ) break;
		//      if ( distance <= 3.5 ) {
		//       RNA_binding_res_in_contact.push_back( i );
		//       protein_res_contacted = true;
		//       break; // only need 1 for it to count as in contact
		//      }
		//     }
		//     if ( protein_res_contacted ) break;
		//    }
		//    if ( protein_res_contacted ) break;
		//   }
		//  }
		//
		//  core::Real num_contacts = RNA_binding_res_in_contact.size();
		//  TR << TR.Blue << tag << " contacts " << num_contacts << "/" << option[ RNA_binding_res ]().size() << " predicted RNA binding residues" << std::endl;

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	count_contacts();

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
		NEW_OPT( RNA_binding_res, "residues that are predicted to bind RNA", 1);
		NEW_OPT( native_RNA_binding_res, "get the native RNA binding residues", false);

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
		return -1;
	}
}
