// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/PackDaemon.cxxtest.hh
/// @brief  test suite for PackDaemon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/pack_daemon/PackDaemon.hh>
#include <protocols/pack_daemon/EntityCorrespondence.hh>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Core headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>
#include <core/pack/interaction_graph/DoubleDensePDInteractionGraph.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSubsets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocols headers
#include <protocols/genetic_algorithm/Entity.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>


//static basic::Tracer TR("PackDaemonTest.cxxtest");

using namespace core;
using namespace protocols::pack_daemon;
using namespace protocols::genetic_algorithm;
using namespace core::chemical;
using namespace core::scoring;
using namespace core::pack;
using namespace core::pack::interaction_graph;
using namespace core::pack::rotamer_set;


class PackDaemonTest : public CxxTest::TestSuite
{
public:
	typedef core::pose::PoseOP PoseOP;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::TaskFactory TaskFactory;
	typedef protocols::pack_daemon::PackDaemon PackDaemon;

public:
	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit" );
	}

	void initialize_default_pack_daemon( PackDaemon & daemon ) {

		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		daemon.set_pose_and_task( *trpcage, *task );
		ScoreFunctionOP sfxn = get_score_function();
		daemon.set_score_function( *sfxn );

		EntityCorrespondenceOP ec( new EntityCorrespondence );
		ec->set_pose( trpcage );
		ec->set_num_entities( 2 );
		ec->add_resid_to_entity_list( 1, 12 );
		ec->add_resid_to_entity_list( 2, 13 );

		daemon.set_entity_correspondence( *ec );

		daemon.setup();

		TS_ASSERT( daemon.correspondence()->entity_for_residue( 11 ) == 0 );
		TS_ASSERT( daemon.correspondence()->entity_for_residue( 12 ) == 1 );
		TS_ASSERT( daemon.correspondence()->entity_for_residue( 13 ) == 2 );
	}

	void test_pack_daemon_ctor() {
		using namespace core;
		using namespace protocols::pack_daemon;
		using namespace protocols::genetic_algorithm;
		using namespace core::chemical;
		using namespace core::scoring;
		using namespace core::pack::rotamer_set;

		PackDaemon daemon;
		initialize_default_pack_daemon( daemon );

		Entity ent( "traits AA:1:F AA:2:R fitness 0.0" );

		RotamerSetsCOP rotsets = daemon.rot_sets();
		/*for ( Size ii = 1; ii <= 3; ++ii ) {
		Size ii_resid = rotsets->moltenres_2_resid( ii );
		Size ii_rot_offset = rotsets->nrotamer_offset_for_moltenres( ii );
		std::cout << "Residue " << ii_resid << " rotamers:";
		RotamerSetCOP ii_rotset = rotsets->rotamer_set_for_moltenresidue( ii );
		for ( Size jj = 1; jj <= ii_rotset->num_rotamers(); ++jj ) {
		std::cout << " (" << ii_rot_offset + jj << "," << ii_rotset->rotamer( jj )->aa() << ")";
		}
		std::cout << std::endl;
		}*/

		utility::vector0< Size > rots_to_pack = daemon.select_rotamer_subset( ent );
		for ( Size ii = 0; ii < rots_to_pack.size(); ++ii ) {
			Size ii_rot = rots_to_pack[ ii ];
			Size ii_moltenres = rotsets->moltenres_for_rotamer( ii_rot );
			Size ii_resid     = rotsets->moltenres_2_resid( ii_moltenres );
			Size ii_entity_id = daemon.correspondence()->entity_for_residue( ii_resid );
			Size ii_local_rot = rotsets->rotid_on_moltenresidue( ii_rot );
			//std::cout << "Rot to pack[" << ii << "] = " << rots_to_pack[ ii ] << " " << ii_moltenres << " " << ii_local_rot << std::endl;
			if ( ii_entity_id == 0 ) continue;
			if ( ii_entity_id == 1 ) {
				TS_ASSERT( rotsets->rotamer_set_for_moltenresidue( ii_moltenres )->rotamer( ii_local_rot )->aa() == aa_phe );
			} else if ( ii_entity_id == 2 ) {
				TS_ASSERT( rotsets->rotamer_set_for_moltenresidue( ii_moltenres )->rotamer( ii_local_rot )->aa() == aa_arg );
			}
		}
		//std::cout << std::endl;
	}

	void test_pack_daemon_repack() {
		using namespace protocols::genetic_algorithm;

		PackDaemon daemon;
		initialize_default_pack_daemon( daemon );
		Entity ent( "traits AA:1:F AA:2:R fitness 0.0" );

		daemon.compute_energy_for_assignment( ent );
		PackDaemon::RotamerAssignmentAndEnergy const & assignment = daemon.last_assignment();
		//std::cout << "Total energy: " << assignment.second << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) {
		// std::cout << "Moltenres " << ii << " with rotamer " << assignment.first[ ii ] << std::endl;
		//}
		//Moltenres 1 with rotamer 1
		//Moltenres 2 with rotamer 5
		//Moltenres 3 with rotamer 42

		//std::cout << "Total energy: " << assignment.second << std::endl;
		for ( Size ii = 1; ii <= 3; ++ii ) {
			//std::cout << "Moltenres " << ii << " with rotamer " << assignment.first[ ii ] << std::endl;
			switch ( ii ) {
			case 1 : TS_ASSERT( assignment.first[ ii ] == 1 ); break;
			case 2 : TS_ASSERT( assignment.first[ ii ] == 42 ); break;
			case 3 : TS_ASSERT( assignment.first[ ii ] == 90 ); break;
			}
		}

	}

	void test_daemon_set_setup() {
		using namespace protocols::genetic_algorithm;
		DaemonSet ds;
		ScoreFunctionOP sfxn = get_score_function();
		ds.set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::string entity_corr_string( "1 12 A\n2 6 A\n1 18 A\n" );
		std::string secondary_resfile_string( "NATRO\nstart\n9 A NATAA\n13 A NATAA" );

		std::istringstream corr_resfile( corr_resfile_string );
		std::istringstream entity_corr( entity_corr_string );
		std::istringstream sec_resfile( secondary_resfile_string );
		PoseOP trpcage = create_trpcage_ideal_poseop();

		ds.set_entity_resfile( corr_resfile, "unnamed" );
		ds.add_pack_daemon( 1, "1l2y.pdb", *trpcage, "entity_corr_string", entity_corr, "secondary_resfile_string",sec_resfile );
		ds.setup_daemons();

		DaemonSet::ConstDaemonList daemons = ds.daemons();
		/*for ( std::list< std::pair< core::Size, PackDaemonCOP > >::const_iterator
		iter = daemons.begin(), iter_end = daemons.end(); iter != iter_end; ++iter ) {
		PackDaemon const & daemon = *( iter->second );
		core::pack::rotamer_set::RotamerSetsCOP rotsets = daemon.rot_sets();
		for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		for ( Size jj = 1; jj <= rotsets->rotamer_set_for_moltenresidue( ii )->num_rotamers(); ++jj ) {
		std::cout << "Rotamer on residue " << ii << " #" << jj << " " << rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( jj )->aa() << std::endl;

		}
		}
		}*/

		Entity ent( "traits AA:1:P AA:2:F fitness 0.0" );

		DaemonSet::StateEsAndNPDs energies = ds.compute_energy_for_assignment( ent );
		for ( std::list< std::pair< core::Size, core::Real > >::const_iterator
				iter = energies.first.begin(), iter_end = energies.first.end();
				iter != iter_end; ++iter ) {
			//std::cout << "Energies: " << iter->first << " " << iter->second << std::endl;
		}

		//std::list< std::pair< core::Size, PackDaemonCOP > > daemons = ds.daemons();
		for ( DaemonSet::ConstDaemonListIter iter = daemons.begin(),
				iter_end = daemons.end(); iter != iter_end; ++iter ) {
			PackDaemon const & daemon = *( iter->second );

			PackDaemon::RotamerAssignmentAndEnergy const & assignment = daemon.last_assignment();
			//std::cout << "Total energy: " << assignment.second << std::endl;
			for ( Size ii = 1; ii <= 3; ++ii ) {
				//std::cout << "Moltenres " << ii << " with rotamer " << assignment.first[ ii ] << std::endl;
				switch ( ii ) {
				case 1 : TS_ASSERT( assignment.first[ ii ] == 22  ); break;
				case 2 : TS_ASSERT( assignment.first[ ii ] == 277 ); break;
				case 3 : TS_ASSERT( assignment.first[ ii ] == 285 ); break;
				}
			}
		}

		Entity bad_ent( "traits AA:1:P AA:2:R fitness 0.0" ); // There are no arginine rotamers for res 2.  This entity will throw an exception
		try {
			ds.compute_energy_for_assignment( bad_ent );
			TS_ASSERT( false ); /// This should not be reached.
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			//std::cout << excn.msg() << std::endl;
			TS_ASSERT( excn.msg() == "Failed to find any rotamers for residue 6 corresponding to entity 2 when looking for ARG rotamers." );
		}

	}

	void test_DenseIGRepacker_create_dense_ig() {
		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_phe ] = allowed_aas[ aa_leu ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 7 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
		ScoreFunctionOP sfxn = get_score_function();

		RotamerSetsOP rot_sets( new RotamerSets );
		InteractionGraphBaseOP ig;
		core::pack::pack_rotamers_setup( *trpcage, *sfxn, task, rot_sets, ig );

		PrecomputedPairEnergiesInteractionGraphOP precomp_ig
			= utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph > ( ig );

		TS_ASSERT( precomp_ig != 0 );
		if ( ! precomp_ig ) return;

		utility::vector0< int > rtp; // Rotamers To Pack.
		for ( Size ii = 1; ii <= 3; ++ii ) {
			for ( Size jj = 1; jj <= rot_sets->nrotamers_for_moltenres( ii ); ++jj ) {
				if ( ii == 1 ) {
					// look for phe rots
					if ( rot_sets->rotamer_for_moltenres( ii, jj )->aa() == aa_leu ) {
						rtp.push_back( rot_sets->nrotamer_offset_for_moltenres( ii ) + jj );
					}
				} else if ( ii == 2 ) {
					// look for leu rots
					if ( rot_sets->rotamer_for_moltenres( ii, jj )->aa() == aa_phe ) {
						rtp.push_back( rot_sets->nrotamer_offset_for_moltenres( ii ) + jj );
					}
				} else {
					// look for arg rots
					if ( rot_sets->rotamer_for_moltenres( ii, jj )->aa() == aa_arg ) {
						rtp.push_back( rot_sets->nrotamer_offset_for_moltenres( ii ) + jj );
					}
				}
			}
		}

		bool all_good = true;
		/* Uncomment to regenerate the rotamer mapping.
		for ( Size ii = 0; ii < rtp.size(); ++ii ) {
		std::cout << "TS_ASSERT( rtp[ " << ii << " ] == " << rtp[ ii ] << ");";
		std::cout << " if ( rtp[ " << ii << " ] != " << rtp[ ii ] << ") all_good = false;" << std::endl;
		}*/

		TS_ASSERT( rtp[ 0 ] == 2); if ( rtp[ 0 ] != 2 ) all_good = false;
		TS_ASSERT( rtp[ 1 ] == 3); if ( rtp[ 1 ] != 3 ) all_good = false;
		TS_ASSERT( rtp[ 2 ] == 25); if ( rtp[ 2 ] != 25 ) all_good = false;
		TS_ASSERT( rtp[ 3 ] == 26); if ( rtp[ 3 ] != 26 ) all_good = false;
		TS_ASSERT( rtp[ 4 ] == 27); if ( rtp[ 4 ] != 27 ) all_good = false;
		TS_ASSERT( rtp[ 5 ] == 28); if ( rtp[ 5 ] != 28 ) all_good = false;
		TS_ASSERT( rtp[ 6 ] == 66); if ( rtp[ 6 ] != 66 ) all_good = false;
		TS_ASSERT( rtp[ 7 ] == 67); if ( rtp[ 7 ] != 67 ) all_good = false;
		TS_ASSERT( rtp[ 8 ] == 68); if ( rtp[ 8 ] != 68 ) all_good = false;
		TS_ASSERT( rtp[ 9 ] == 69); if ( rtp[ 9 ] != 69 ) all_good = false;
		TS_ASSERT( rtp[ 10 ] == 70); if ( rtp[ 10 ] != 70 ) all_good = false;
		TS_ASSERT( rtp[ 11 ] == 71); if ( rtp[ 11 ] != 71 ) all_good = false;
		TS_ASSERT( rtp[ 12 ] == 72); if ( rtp[ 12 ] != 72 ) all_good = false;
		TS_ASSERT( rtp[ 13 ] == 73); if ( rtp[ 13 ] != 73 ) all_good = false;
		TS_ASSERT( rtp[ 14 ] == 74); if ( rtp[ 14 ] != 74 ) all_good = false;
		TS_ASSERT( rtp[ 15 ] == 75); if ( rtp[ 15 ] != 75 ) all_good = false;
		TS_ASSERT( rtp[ 16 ] == 76); if ( rtp[ 16 ] != 76 ) all_good = false;
		TS_ASSERT( rtp[ 17 ] == 77); if ( rtp[ 17 ] != 77 ) all_good = false;
		TS_ASSERT( rtp[ 18 ] == 78); if ( rtp[ 18 ] != 78 ) all_good = false;
		TS_ASSERT( rtp[ 19 ] == 79); if ( rtp[ 19 ] != 79 ) all_good = false;
		TS_ASSERT( rtp[ 20 ] == 80); if ( rtp[ 20 ] != 80 ) all_good = false;
		TS_ASSERT( rtp[ 21 ] == 81); if ( rtp[ 21 ] != 81 ) all_good = false;
		TS_ASSERT( rtp[ 22 ] == 82); if ( rtp[ 22 ] != 82 ) all_good = false;

		TS_ASSERT( all_good );
		if ( ! all_good ) return;

		DenseIGRepackerOP repacker( new DenseIGRepacker( trpcage, task, precomp_ig, rot_sets ) );
		RotamerSubsetsOP rotsubset = repacker->create_rotamer_subsets_from_rot_to_pack( rtp );
		DensePDInteractionGraphOP denseig = repacker->create_dense_pdig_from_rot_to_pack( rtp, rotsubset );

		ObjexxFCL::FArray1D_int network_state_orig( 3 );
		ObjexxFCL::FArray1D_int network_state_sset( 3 );

		precomp_ig->prepare_for_simulated_annealing();
		denseig->prepare_for_simulated_annealing();

		network_state_sset( 1 ) = 1; network_state_sset( 2 ) = 4; network_state_sset( 3 ) = 18;
		for ( Size ii = 1; ii <= 3; ++ii ) {
			network_state_orig( ii ) = rot_sets->rotid_on_moltenresidue( rtp[ network_state_sset( ii ) - 1 ] );
			network_state_sset( ii ) = rotsubset->rotid_on_moltenresidue( network_state_sset( ii ) );
		}

		//std::cout << "orig: " << network_state_orig( 1 ) << " " << network_state_orig( 2 ) << " " << network_state_orig( 3 ) << std::endl;
		//std::cout << "sset: " << network_state_sset( 1 ) << " " << network_state_sset( 2 ) << " " << network_state_sset( 3 ) << std::endl;

		TS_ASSERT_DELTA( denseig->set_network_state( network_state_sset ), precomp_ig->set_network_state( network_state_orig ), 1e-5 );

		network_state_sset( 1 ) = 2; network_state_sset( 2 ) = 5; network_state_sset( 3 ) = 7;
		for ( Size ii = 1; ii <= 3; ++ii ) {
			network_state_orig( ii ) = rot_sets->rotid_on_moltenresidue( rtp[ network_state_sset( ii ) - 1 ] );
			network_state_sset( ii ) = rotsubset->rotid_on_moltenresidue( network_state_sset( ii ) );
		}

		//std::cout << "orig: " << network_state_orig( 1 ) << " " << network_state_orig( 2 ) << " " << network_state_orig( 3 ) << std::endl;
		//std::cout << "sset: " << network_state_sset( 1 ) << " " << network_state_sset( 2 ) << " " << network_state_sset( 3 ) << std::endl;

		TS_ASSERT_DELTA( denseig->set_network_state( network_state_sset ), precomp_ig->set_network_state( network_state_orig ), 1e-5 );


	}

	void test_DoubleDenseIGRepacker_create_double_dense_ig() {
		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_phe ] = allowed_aas[ aa_leu ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 7 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
		ScoreFunctionOP sfxn = get_score_function();

		RotamerSetsOP rot_sets( new RotamerSets );
		InteractionGraphBaseOP ig;
		core::pack::pack_rotamers_setup( *trpcage, *sfxn, task, rot_sets, ig );

		PrecomputedPairEnergiesInteractionGraphOP precomp_ig
			= utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph > ( ig );

		TS_ASSERT( precomp_ig != 0 );
		if ( ! precomp_ig ) return;

		utility::vector0< int > rtp; // Rotamers To Pack.
		for ( Size ii = 1; ii <= 3; ++ii ) {
			for ( Size jj = 1; jj <= rot_sets->nrotamers_for_moltenres( ii ); ++jj ) {
				if ( ii == 1 ) {
					// look for phe rots
					if ( rot_sets->rotamer_for_moltenres( ii, jj )->aa() == aa_leu ) {
						rtp.push_back( rot_sets->nrotamer_offset_for_moltenres( ii ) + jj );
					}
				} else if ( ii == 2 ) {
					// look for leu rots
					if ( rot_sets->rotamer_for_moltenres( ii, jj )->aa() == aa_phe ) {
						rtp.push_back( rot_sets->nrotamer_offset_for_moltenres( ii ) + jj );
					}
				} else {
					// look for arg rots
					if ( rot_sets->rotamer_for_moltenres( ii, jj )->aa() == aa_arg ) {
						rtp.push_back( rot_sets->nrotamer_offset_for_moltenres( ii ) + jj );
					}
				}
			}
		}

		bool all_good = true;
		/* Uncomment to regenerate the rotamer mapping.
		for ( Size ii = 0; ii < rtp.size(); ++ii ) {
		std::cout << "TS_ASSERT( rtp[ " << ii << " ] == " << rtp[ ii ] << ");";
		std::cout << " if ( rtp[ " << ii << " ] != " << rtp[ ii ] << ") all_good = false;" << std::endl;
		}*/

		TS_ASSERT( rtp[ 0 ] == 2); if ( rtp[ 0 ] != 2 ) all_good = false;
		TS_ASSERT( rtp[ 1 ] == 3); if ( rtp[ 1 ] != 3 ) all_good = false;
		TS_ASSERT( rtp[ 2 ] == 25); if ( rtp[ 2 ] != 25 ) all_good = false;
		TS_ASSERT( rtp[ 3 ] == 26); if ( rtp[ 3 ] != 26 ) all_good = false;
		TS_ASSERT( rtp[ 4 ] == 27); if ( rtp[ 4 ] != 27 ) all_good = false;
		TS_ASSERT( rtp[ 5 ] == 28); if ( rtp[ 5 ] != 28 ) all_good = false;
		TS_ASSERT( rtp[ 6 ] == 66); if ( rtp[ 6 ] != 66 ) all_good = false;
		TS_ASSERT( rtp[ 7 ] == 67); if ( rtp[ 7 ] != 67 ) all_good = false;
		TS_ASSERT( rtp[ 8 ] == 68); if ( rtp[ 8 ] != 68 ) all_good = false;
		TS_ASSERT( rtp[ 9 ] == 69); if ( rtp[ 9 ] != 69 ) all_good = false;
		TS_ASSERT( rtp[ 10 ] == 70); if ( rtp[ 10 ] != 70 ) all_good = false;
		TS_ASSERT( rtp[ 11 ] == 71); if ( rtp[ 11 ] != 71 ) all_good = false;
		TS_ASSERT( rtp[ 12 ] == 72); if ( rtp[ 12 ] != 72 ) all_good = false;
		TS_ASSERT( rtp[ 13 ] == 73); if ( rtp[ 13 ] != 73 ) all_good = false;
		TS_ASSERT( rtp[ 14 ] == 74); if ( rtp[ 14 ] != 74 ) all_good = false;
		TS_ASSERT( rtp[ 15 ] == 75); if ( rtp[ 15 ] != 75 ) all_good = false;
		TS_ASSERT( rtp[ 16 ] == 76); if ( rtp[ 16 ] != 76 ) all_good = false;
		TS_ASSERT( rtp[ 17 ] == 77); if ( rtp[ 17 ] != 77 ) all_good = false;
		TS_ASSERT( rtp[ 18 ] == 78); if ( rtp[ 18 ] != 78 ) all_good = false;
		TS_ASSERT( rtp[ 19 ] == 79); if ( rtp[ 19 ] != 79 ) all_good = false;
		TS_ASSERT( rtp[ 20 ] == 80); if ( rtp[ 20 ] != 80 ) all_good = false;
		TS_ASSERT( rtp[ 21 ] == 81); if ( rtp[ 21 ] != 81 ) all_good = false;
		TS_ASSERT( rtp[ 22 ] == 82); if ( rtp[ 22 ] != 82 ) all_good = false;

		TS_ASSERT( all_good );
		if ( ! all_good ) return;

		DoubleDenseIGRepackerOP repacker( new DoubleDenseIGRepacker( trpcage, task, precomp_ig, rot_sets ) );
		RotamerSubsetsOP rotsubset = repacker->create_rotamer_subsets_from_rot_to_pack( rtp );
		DoubleDensePDInteractionGraphOP double_dense_ig = repacker->create_dense_pdig_from_rot_to_pack( rtp, rotsubset );

		ObjexxFCL::FArray1D_int network_state_orig( 3 );
		ObjexxFCL::FArray1D_int network_state_sset( 3 );

		precomp_ig->prepare_for_simulated_annealing();
		double_dense_ig->prepare_for_simulated_annealing();

		precomp_ig->blanket_assign_state_0();
		double_dense_ig->blanket_assign_state_0();

		network_state_sset( 1 ) = 1; network_state_sset( 2 ) = 4; network_state_sset( 3 ) = 18;
		for ( Size ii = 1; ii <= 3; ++ii ) {
			network_state_orig( ii ) = rot_sets->rotid_on_moltenresidue( rtp[ network_state_sset( ii ) - 1 ] );
			network_state_sset( ii ) = rotsubset->rotid_on_moltenresidue( network_state_sset( ii ) );
		}

		//std::cout << "orig: " << network_state_orig( 1 ) << " " << network_state_orig( 2 ) << " " << network_state_orig( 3 ) << std::endl;
		//std::cout << "sset: " << network_state_sset( 1 ) << " " << network_state_sset( 2 ) << " " << network_state_sset( 3 ) << std::endl;

		TS_ASSERT_DELTA( double_dense_ig->set_network_state( network_state_sset ), precomp_ig->set_network_state( network_state_orig ), 1e-5 );

		PackerEnergy dd_deltaE( 0.0 ), sparse_deltaE( 0.0 ), dd_prevres_totE( 0.0 ), sparse_prevres_totE( 0.0 );

		double_dense_ig->consider_substitution( 3, rotsubset->rotid_on_moltenresidue( 6 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 3, rot_sets->rotid_on_moltenresidue( rtp[ 6 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

		double_dense_ig->consider_substitution( 3, rotsubset->rotid_on_moltenresidue( 7 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 3, rot_sets->rotid_on_moltenresidue( rtp[ 7 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

		double_dense_ig->consider_substitution( 3, rotsubset->rotid_on_moltenresidue( 8 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 3, rot_sets->rotid_on_moltenresidue( rtp[ 8 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

		double_dense_ig->consider_substitution( 3, rotsubset->rotid_on_moltenresidue( 9 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 3, rot_sets->rotid_on_moltenresidue( rtp[ 9 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

		double_dense_ig->consider_substitution( 3, rotsubset->rotid_on_moltenresidue( 10 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 3, rot_sets->rotid_on_moltenresidue( rtp[ 10 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );


		network_state_sset( 1 ) = 2; network_state_sset( 2 ) = 5; network_state_sset( 3 ) = 7;
		for ( Size ii = 1; ii <= 3; ++ii ) {
			network_state_orig( ii ) = rot_sets->rotid_on_moltenresidue( rtp[ network_state_sset( ii ) - 1 ] );
			network_state_sset( ii ) = rotsubset->rotid_on_moltenresidue( network_state_sset( ii ) );
		}

		//std::cout << "orig: " << network_state_orig( 1 ) << " " << network_state_orig( 2 ) << " " << network_state_orig( 3 ) << std::endl;
		//std::cout << "sset: " << network_state_sset( 1 ) << " " << network_state_sset( 2 ) << " " << network_state_sset( 3 ) << std::endl;

		TS_ASSERT_DELTA( double_dense_ig->set_network_state( network_state_sset ), precomp_ig->set_network_state( network_state_orig ), 1e-5 );

		double_dense_ig->consider_substitution( 2, rotsubset->rotid_on_moltenresidue( 2 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 2, rot_sets->rotid_on_moltenresidue( rtp[ 2 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

		double_dense_ig->consider_substitution( 2, rotsubset->rotid_on_moltenresidue( 3 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 2, rot_sets->rotid_on_moltenresidue( rtp[ 3 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

		double_dense_ig->consider_substitution( 2, rotsubset->rotid_on_moltenresidue( 4 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 2, rot_sets->rotid_on_moltenresidue( rtp[ 4 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

		double_dense_ig->consider_substitution( 2, rotsubset->rotid_on_moltenresidue( 5 + 1 ), dd_deltaE, dd_prevres_totE );
		precomp_ig->consider_substitution( 2, rot_sets->rotid_on_moltenresidue( rtp[ 5 ] ), sparse_deltaE, sparse_prevres_totE );

		TS_ASSERT_DELTA( dd_deltaE, sparse_deltaE, 1e-6 );
		TS_ASSERT_DELTA( dd_prevres_totE, sparse_prevres_totE, 1e-6 );

	}


};


