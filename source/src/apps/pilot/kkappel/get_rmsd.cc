// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Calculate RMSD with option of align residues and RMSD residues

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
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
using namespace core::import_pose;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

OPT_KEY( Boolean, dump_structures )
OPT_KEY( Boolean, rmsd_nosuper )
OPT_KEY( Boolean, backbone_rmsd )
OPT_KEY( Boolean, heavy_atom_rmsd )
OPT_KEY( Boolean, protein_align )
OPT_KEY( Boolean, rna_rmsd )
OPT_KEY( IntegerVector, rmsd_res )
OPT_KEY( IntegerVector, native_rmsd_res )
OPT_KEY( IntegerVector, align_residues )
OPT_KEY( IntegerVector, native_align_residues )

///////////////////////////////////////////////////////////////////////////////
void
get_rmsd()
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
	using namespace protocols::stepwise::modeler;
	using namespace protocols::stepwise::setup;

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

	// native pose setup
	pose::PoseOP native_pose;
	if ( option[ in::file::native ].user() ) {
		native_pose = get_pdb_with_full_model_info( option[ in::file::native ](), rsd_set );
		native_pose->dump_pdb("native_pose_dump.pdb");
	}


	// score function setup
	core::scoring::ScoreFunctionOP scorefxn;
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
	}

	FullModelInfoOP my_model;
	utility::vector1< pose::PoseOP > other_poses;

	pose::Pose pose,start_pose;

	Size i( 0 );

	// Silent file output setup
	std::string const silent_file = option[ out::file::silent  ]();
	SilentFileOptions opts; // initialized from the command line
	SilentFileData silent_file_data( opts );


	//if ( native_exists ) ( *scorefxn)( *native_pose );

	while ( input->has_another_pose() ) {

		input->fill_pose( pose, *rsd_set );
		i++;

		//if ( !option[ in::file::silent ].user() ) cleanup( pose );

		// tag
		std::string tag = tag_from_pose( pose );

		// if align over protein

		utility::vector1< core::Size > native_super_residues, super_residues;
		if ( !option[ rmsd_nosuper ]() ) {
			if ( option[ protein_align ]() ) {
				for ( core::Size i = 1; i<= pose.total_residue(); ++i ) {
					if ( pose.residue( i).is_protein() ) {
						super_residues.push_back( i );
					}
				}
				for ( core::Size i = 1; i<= native_pose->total_residue(); ++i ) {
					if ( native_pose->residue( i ).is_protein() ) {
						native_super_residues.push_back( i );
					}
				}
			} else {
				native_super_residues = option[ native_align_residues ]();
				super_residues = option[ align_residues ]();
			}
		}

		// if calculate rmsd over RNA:
		utility::vector1< core::Size > native_rmsd_residues, rmsd_residues;
		if ( option[ rna_rmsd ]() ) {
			for ( core::Size i = 1; i<= pose.total_residue(); ++i ) {
				if ( pose.residue( i ).is_RNA() ) {
					rmsd_residues.push_back( i );
				}
			}
			for ( core::Size i = 1; i<= native_pose->total_residue(); ++i ) {
				if ( native_pose->residue( i ).is_RNA() ) {
					native_rmsd_residues.push_back( i );
				}
			}
			//std::cout << "RMSD RESIDUES: " << rmsd_residues << std::endl;
			//std::cout << "NATIVE RMSD RESIDUES: " << native_rmsd_residues << std::endl;
		} else {
			native_rmsd_residues = option[ native_rmsd_res ]();
			rmsd_residues = option[ rmsd_res ]();
		}

		// calculate RMSD
		Real rmsd;
		Real align_rmsd( 0.0 );
		if ( !option[ rmsd_nosuper ]() ) {
			// first do the superposition
			core::id::AtomID_Map< id::AtomID > atom_map;
			core::pose::initialize_atomid_map( atom_map, pose, id::AtomID::BOGUS_ATOM_ID() );
			for ( core::Size i = 1; i <= native_super_residues.size(); ++i ) {
				id::AtomID const id1( pose.residue( super_residues[ i ] ).atom_index("CA"), super_residues[ i ]);
				id::AtomID const id2( native_pose->residue( native_super_residues[ i ]).atom_index("CA"), native_super_residues[ i ] );
				atom_map[ id1 ] = id2;
			}

			align_rmsd = superimpose_pose( pose, *native_pose, atom_map );
		}

		// RNA backbone atoms:
		utility::vector1< std::string > RNA_backbone_atoms;
		RNA_backbone_atoms.push_back(" P  ");
		RNA_backbone_atoms.push_back(" OP1");
		RNA_backbone_atoms.push_back(" OP2");
		RNA_backbone_atoms.push_back(" O5'");
		RNA_backbone_atoms.push_back(" C5'");
		RNA_backbone_atoms.push_back(" C4'");
		RNA_backbone_atoms.push_back(" O4'");
		RNA_backbone_atoms.push_back(" C3'");
		RNA_backbone_atoms.push_back(" O3'");
		RNA_backbone_atoms.push_back(" C1'");
		RNA_backbone_atoms.push_back(" C2'");
		RNA_backbone_atoms.push_back(" O2'");

		// do backbone RMSD
		// then get the rmsd
		//core::id::AtomID_Map< id::AtomID > atom_map_rms;
		std::map<core::id::AtomID, core::id::AtomID> atom_map_rms;
		//core::pose::initialize_atomid_map( atom_map_rms, pose, id::BOGUS_ATOM_ID );
		if ( option[ backbone_rmsd ]() ) {
			for ( core::Size i = 1; i<= native_rmsd_residues.size(); ++i ) {
				for ( core::Size j = 1; j<= RNA_backbone_atoms.size(); ++j ) {
					std::string name = RNA_backbone_atoms[ j ];
					//std::cout << "CHECKING native residue " << native_rmsd_residues[i] << " rmsd residue " << rmsd_residues[i] << " name " << name <<std::endl;
					if ( !pose.residue(rmsd_residues[i]).has( name ) ) continue;
					if ( !native_pose->residue(native_rmsd_residues[i]).has( name ) ) continue;
					id::AtomID const id1( pose.residue( rmsd_residues[ i ] ).atom_index(name), rmsd_residues[ i ]);
					id::AtomID const id2( native_pose->residue( native_rmsd_residues[ i ]).atom_index(name), native_rmsd_residues[ i ] );
					atom_map_rms[ id1 ] = id2;
					//std::cout << "added native residue " << native_rmsd_residues[i] << " rmsd residue " << rmsd_residues[i] << " name " << name <<std::endl;
				}
				//id::AtomID const id1( pose.residue( rmsd_residues[ i ] ).atom_index(" P  "), rmsd_residues[ i ]);
				//id::AtomID const id2( native_pose->residue( native_rmsd_residues[ i ]).atom_index(" P  "), native_rmsd_residues[ i ] );
				//atom_map_rms[ id1 ] = id2;
			}
		} else if ( option[ heavy_atom_rmsd ] () ) {
			// do heavy atom RMSD
			for ( core::Size i = 1; i<= native_rmsd_residues.size(); ++i ) {
				core::conformation::Residue const & rsd1 = native_pose->residue( native_rmsd_residues[i] );
				core::conformation::Residue const & rsd2 = pose.residue( rmsd_residues[i] );

				for ( Size j = 1; j<=rsd1.nheavyatoms(); ++j ) {
					std::string name( rsd1.atom_name( j ) );
					//std::cout << "CHECKING native residue " << native_rmsd_residues[i] << " rmsd residue " << rmsd_residues[i] << " name " << name <<std::endl;
					if ( !rsd2.has( name ) ) continue;
					if ( rsd1.is_virtual(j) ) continue;
					Size const j2( rsd2.atom_index( name ) );
					if ( rsd2.is_virtual(j2) ) continue;
					id::AtomID const id1( rsd1.atom_index( name ), native_rmsd_residues[ i ] );
					id::AtomID const id2( rsd2.atom_index( name ), rmsd_residues[ i ] );
					atom_map_rms[ id2 ] = id1;
					//atom_map_rms[ id1 ] = id2;
					//std::cout << "added native residue " << native_rmsd_residues[i] << " rmsd residue " << rmsd_residues[i] << " name " << name <<std::endl;
				}
			}
		}

		rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose, atom_map_rms );
		//rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose, atom_map_rms, rmsd_residues );

		//  if ( native_exists ) {
		//   Real rmsd;
		//   if ( option[ rmsd_nosuper ]() ) {
		//    if ( option[ rmsd_residues ].user() ) {
		//     rmsd = protocols::stepwise::modeler::align::get_rmsd( pose, *native_pose, option[ rmsd_residues ]() );
		//    } else {
		//     rmsd      = all_atom_rmsd_nosuper( *native_pose, pose );
		//    }
		//   } else {
		//    //Real const rmsd      = all_atom_rmsd( *native_pose, pose );
		//    rmsd = protocols::stepwise::modeler::align::superimpose_with_stepwise_aligner( pose, *native_pose, option[ OptionKeys::stepwise::superimpose_over_all ]() );
		//    //Real const rmsd = protocols::stepwise::modeler::align::superimpose_with_stepwise_aligner( pose, *native_pose, option[ OptionKeys::stepwise::superimpose_over_all ]() );
		//   }
		//   std::cout << "All atom rmsd: " << tag << " " << rmsd << std::endl;
		//
		//  }
		std::cout << "Align rmsd: " << tag << " " << align_rmsd << std::endl;
		std::cout << "All atom rmsd: " << tag << " " << rmsd << std::endl;

		// write it to the silent file
		BinarySilentStruct s( opts, pose, tag );
		// could also score it here?
		s.add_energy( "rms", rmsd );
		s.add_energy( "align_rms", align_rmsd );
		silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );

		if ( option[ dump_structures ]() && tag.find( ".pdb" ) != std::string::npos ) {
			std::string out_pdb_file = utility::replace_in( tag, ".pdb", ".sup.pdb" );
			std::cout << "Creating: " << out_pdb_file << std::endl;
			pose.dump_pdb( out_pdb_file );
		}

		//// backbone atoms
		//for ( Size i = 1; i <= native_pose->residue( 2 ).nheavyatoms(); ++i ) {
		// if (native_pose->residue( 2 ).atom_is_hydrogen( i )) continue;
		// std::cout << native_pose->residue(2).atom_name( i ) << " is backbone: " << native_pose->residue(2).atom_is_backbone(2) << std::endl;
		//}

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	get_rmsd();

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
		NEW_OPT( rmsd_nosuper, "Calculate rmsd without superposition first", false);
		NEW_OPT( rna_rmsd, "Calculate rmsd over rna", true);
		NEW_OPT( backbone_rmsd, "Calculate RNA heavy backbone atoms", false);
		NEW_OPT( heavy_atom_rmsd, "Calculate rmsd over all heavy atoms", false);
		NEW_OPT( protein_align, "align structures over protein residues", true);
		NEW_OPT( dump_structures, "dump superimposed structures", false);
		NEW_OPT( rmsd_res, "residues to calculate rmsd for", 1);
		NEW_OPT( native_rmsd_res, "native residues to calculate rmsd for", 1);
		NEW_OPT( align_residues, "residues to align over, default all", 1);
		NEW_OPT( native_align_residues, "native residues to align over, default all", 1);

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
