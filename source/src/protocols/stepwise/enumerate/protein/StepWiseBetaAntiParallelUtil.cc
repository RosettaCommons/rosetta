// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseBetaAntiParallelJumpSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseBetaAntiParallelUtil.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/stepwise/enumerate/general/StepWiseClusterer.hh>
#include <utility/io/ozstream.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/NamedAtomID.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/pose/Pose.hh>

#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <numeric/xyz.functions.hh>

#include <ObjexxFCL/string.functions.hh>

#include <iostream>
#include <string>

//Auto Headers
#include <core/io/silent/SilentFileData.hh>
#include <utility/vector1.hh>



using namespace core;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

	////////////////////////////////////////////////////////
	// Move ths out of stepwise_protein_test.cc!!
	////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	void
	do_set_xyz( pose::Pose const & pose, Size const i, pose::Pose & scratch_pose, Size const i_scratch, core::kinematics::Stub const & stub ) {

		using namespace core::id;
		using namespace core::chemical;

		scratch_pose.set_xyz( NamedAtomID( " CA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " CA ", i) ) ) );
		scratch_pose.set_xyz( NamedAtomID( " N  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " N  ", i) ) ) );
		scratch_pose.set_xyz( NamedAtomID( " C  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " C  ", i) ) ) );
		scratch_pose.set_xyz( NamedAtomID( " O  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " O  ", i) ) ) );

		if ( pose.residue_type(i).has( " H  " ) ){
			scratch_pose.set_xyz( NamedAtomID( " H  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " H  ", i) ) ) );
		} else if ( pose.residue_type(i).has( " H1 " ) ){
			scratch_pose.set_xyz( NamedAtomID( " H  ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " H1 ", i) ) ) );
		}

		if ( pose.residue_type(i).has( "1HA " ) ){
			scratch_pose.set_xyz( NamedAtomID( "1HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( "1HA ", i) ) ) );
		} else if ( pose.residue_type(i).has( " HA " ) ){
			scratch_pose.set_xyz( NamedAtomID( "1HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " HA ", i) ) ) );
		}

		if ( pose.residue_type(i).has( "2HA " ) ){
			scratch_pose.set_xyz( NamedAtomID( "2HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( "2HA ", i) ) ) );
		} else if ( pose.residue_type(i).has( " CB " ) ){
			scratch_pose.set_xyz( NamedAtomID( "2HA ", i_scratch), stub.global2local( pose.xyz( NamedAtomID( " CB ", i) ) ) );
		}


	}


	////////////////////////////////////////////////////////
	// Move ths out of stepwise_protein_test.cc!!
	////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	void
	generate_beta_database_test(){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		using namespace core::scoring::hbonds;
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::conformation;
		using namespace core::id;
		using namespace core::io::silent;

		static ScoreFunctionOP scorefxn = getScoreFunction();
		ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		Pose pose, scratch_pose;

		// little pose will be used to output beta geometries.
		std::string sequence = "GG";
		core::pose::make_pose_from_sequence( scratch_pose,	sequence,	*rsd_set );
		FoldTree f( scratch_pose.fold_tree() );

		f.new_jump(1,2,1);
		f.set_jump_atoms( 1, " CA ", " CA ", true /*bKeepStubInResidue. I think this is still protein-centric*/ );
		f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.
		scratch_pose.fold_tree( f );

	pose::remove_variant_type_from_pose_residue( scratch_pose, LOWER_TERMINUS, 1 );
	pose::remove_variant_type_from_pose_residue( scratch_pose, UPPER_TERMINUS, 2 );

		// main loop
		utility::vector1< std::string > in_files = protocols::stepwise::load_s_and_l();

		SilentFileDataOP silent_file_data = new SilentFileData;
		std::string const silent_file( option[ out::file::silent ]() );

		Size count( 0 );
		for ( Size n = 1; n <= in_files.size(); n++ ) {

			std::cout << "-------- " << in_files[n] << "----------" << std::endl;
			import_pose::pose_from_pdb( pose, *rsd_set, in_files[ n ] );

			// Look for beta pairings (antiparallel for now)
			(*scorefxn)(pose);
			HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
			hbond_options->use_hb_env_dep( false );
			HBondSetOP hbond_set( new hbonds::HBondSet( hbond_options ) );

			fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

			std::map< std::pair< Size, Size >, bool > donor_to_acceptor_bb_hb;

			// In antiparallel HB, there are two hbonds between the residue backbones:
			//   C=O ... HN  and a NH ... O=C

			for (Size i = 1; i <= hbond_set->nhbonds(); i++ ) {
				hbonds::HBond const & hbond( hbond_set->hbond( i ) );

				Size const don_res_num = hbond.don_res();
				Size const acc_res_num = hbond.acc_res();

				if ( hbond.don_hatm_is_protein_backbone() && hbond.acc_atm_is_protein_backbone() ) {
					donor_to_acceptor_bb_hb[ std::make_pair( don_res_num, acc_res_num ) ] = true;
				}

			}

			for ( Size i = 1; i <= pose.total_residue(); i++ ){
				for ( Size j = 1; j <= pose.total_residue(); j++ ){
					if ( donor_to_acceptor_bb_hb[ std::make_pair( i, j ) ] &&
							 donor_to_acceptor_bb_hb[ std::make_pair( j, i ) ] ){
						//Create a coordinate system at residue i,
						Stub stub( pose.xyz( NamedAtomID( " CA ", i ) ),
											 pose.xyz( NamedAtomID( " CA ", i ) ),
											 pose.xyz( NamedAtomID( " N  ", i ) ),
											 pose.xyz( NamedAtomID( " C  ", i ) ) );

						// spit out coordinates.
						do_set_xyz( pose, i, scratch_pose, 1, stub );
						do_set_xyz( pose, j, scratch_pose, 2, stub );

						count++;

						std::string const tag = "S_"+ObjexxFCL::string_of( count );
						(*scorefxn)( scratch_pose );
						BinaryProteinSilentStruct s( scratch_pose, tag );

						silent_file_data->write_silent_struct( s, silent_file, false /*write score only*/ );
						silent_file_data->add_structure( s );

					}
				}
			}

		}

		/////////////////////////////////////////////////////////////////////
		// CLUSTER!
		/////////////////////////////////////////////////////////////////////
		protocols::stepwise::enumerate::general::StepWiseClusterer stepwise_clusterer( silent_file_data );
		Size max_decoys( 400 );
		if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
		stepwise_clusterer.set_max_decoys( max_decoys );
		//		stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] ); // false by default
		stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
		stepwise_clusterer.set_score_diff_cut( 10.0 );

		// Do not superimpose inside the clusterer -- trust our superposition over residue 1.
		utility::vector1< Size > calc_rms_residues;
		calc_rms_residues.push_back( 1 );
		calc_rms_residues.push_back( 2 );
		stepwise_clusterer.set_calc_rms_res( calc_rms_residues );

		Real cluster_radius( 0.1 );
		if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
		stepwise_clusterer.set_cluster_radius( cluster_radius	);

		stepwise_clusterer.cluster();

		std::string const silent_file_cluster   = get_file_name( silent_file, "_cluster" );
		stepwise_clusterer.output_silent_file( silent_file_cluster );

		/////////////////////////////////////////////////////////////////////
		// Output jump information
		/////////////////////////////////////////////////////////////////////
		std::string outfile  = option[ out::file::o ];
		utility::io::ozstream out( outfile );

		using namespace protocols::stepwise;
		PoseList const & pose_list = stepwise_clusterer.clustered_pose_list();

		for ( PoseList::const_iterator it = pose_list.begin(); it != pose_list.end(); it++ ){
			pose = *(it->second);
			Stub stub1_forward( pose.xyz( NamedAtomID( " N  ", 1) ),pose.xyz( NamedAtomID( " CA ", 1) ),pose.xyz( NamedAtomID( " C  ", 1) ));
			Stub stub1_reverse( pose.xyz( NamedAtomID( " C  ", 1) ),pose.xyz( NamedAtomID( " CA ", 1) ),pose.xyz( NamedAtomID( " N  ", 1) ));

			Stub stub2_forward( pose.xyz( NamedAtomID( " N  ", 2) ),pose.xyz( NamedAtomID( " CA ", 2) ),pose.xyz( NamedAtomID( " C  ", 2) ));
			Stub stub2_reverse( pose.xyz( NamedAtomID( " C  ", 2) ),pose.xyz( NamedAtomID( " CA ", 2) ),pose.xyz( NamedAtomID( " N  ", 2) ));

			out << "A N N " << 	Jump( stub1_forward, stub2_forward ) << std::endl;
			out << "A N C " << 	Jump( stub1_forward, stub2_reverse ) << std::endl;
			out << "A C N " << 	Jump( stub1_reverse, stub2_forward ) << std::endl;
			out << "A C C " << 	Jump( stub1_reverse, stub2_reverse ) << std::endl;

		}

		out.close();

		std::cout << "Put JUMP transforms in " << outfile << std::endl;

	}



} //protein
} //enumerate
} //stepwise
} //protocols
