// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerTrials.cxxtest.hh
/// @brief  test suite for rotamer_trials
/// @author Florian Richter (floric@u.washington.edu)


// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <utility/graph/Graph.hh>

#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/packer_neighbors.hh>

#include <test/core/init_util.hh>
//#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/types.hh>

#include <basic/basic.hh>

#include <numeric/constants.hh>

#include <test/UTracer.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.rotamer_set.RotamerSet.cxxtest");

using namespace core;

class RotamerSetsTests : public CxxTest::TestSuite
{
	//chemical::ResidueTypeSetCAP residue_set;

public:
	RotamerSetsTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH" );

		//residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief test for rotamer_trials

	void test_rotamer_sets()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack::rotamer_set;
		using namespace pose;

		typedef utility::vector1< core::conformation::ResidueCOP > ResidueCOPs;

		//core::Real const rad_per_deg = numeric::constants::f::degrees_to_radians;
		//core::Real const pi = numeric::constants::f::pi;

		TR << "Beginning RotamerSets test... " << std::endl;

		// init/reset seeds in all RG objects we have to do this inside the test it self function since
		// user could request to run just one single test.
		core::init::init_random_generators(1101, "mt19937");

		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		// read in pose
		Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/1l2y_renameH.pdb" , core::import_pose::PDB_file);

		//pose.dump_pdb("/Users/flo/rosetta/rosetta_source/test/1l2y_rename_H.pdb");

		//score the pose just for the heck of it
		( *scorefxn )( pose );

		// create paker task for rotamer trials
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		//diversify the residue level tasks a bit
		task->nonconst_residue_task( 2 ).restrict_to_repacking();
		task->nonconst_residue_task( 4 ).restrict_to_repacking();
		task->nonconst_residue_task( 5 ).restrict_to_repacking();
		task->nonconst_residue_task( 12 ).restrict_to_repacking();
		task->nonconst_residue_task( 16 ).restrict_to_repacking();
		task->nonconst_residue_task( 17 ).restrict_to_repacking();
		task->nonconst_residue_task( 20 ).restrict_to_repacking();

		task->nonconst_residue_task( 8 ).prevent_repacking();
		task->nonconst_residue_task( 11 ).prevent_repacking();

		task->nonconst_residue_task( 6 ).or_include_current( true );
		task->nonconst_residue_task( 3 ).or_include_current( true );

		task->nonconst_residue_task( 16 ).or_ex1( true );
		task->nonconst_residue_task( 16 ).or_ex2( true );
		task->nonconst_residue_task( 16 ).or_ex3( true );
		task->nonconst_residue_task( 16 ).or_ex4( true );

		task->nonconst_residue_task( 10 ).or_ex1( true );
		task->nonconst_residue_task( 10 ).or_ex2( true );
		task->nonconst_residue_task( 18 ).or_ex1( true );
		task->nonconst_residue_task( 18 ).or_ex2( true );

		task->nonconst_residue_task( 3 ).or_ex1aro( true );
		task->nonconst_residue_task( 3 ).or_ex2aro( true );

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, true );
		keep_aas[ core::chemical::aa_cys ] = false;
		keep_aas[ core::chemical::aa_arg ] = false;
		keep_aas[ core::chemical::aa_glu ] = false;
		keep_aas[ core::chemical::aa_his ] = false;
		keep_aas[ core::chemical::aa_asn ] = false;
		keep_aas[ core::chemical::aa_ile ] = false;
		task->nonconst_residue_task( 18 ).restrict_absent_canonical_aas( keep_aas );

		//aight, let's get a rotamer set corresponding to this task

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, task );
		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );

		//rotsets->dump_pdb( pose, "test_rotamer_sets.pdb");

		TS_ASSERT_EQUALS( rotsets->nrotamers(), /*2314*/2249 );

		//now we'll load the template rotamer set and make sure that all rotamers are equal.
		//we rely on the multimodel pdb reader for file processing

		utility::vector1< core::Size > num_rotamers;
		utility::vector1< core::Size > ref_num_rotamers;

		utility::vector1< ResidueCOPs > ref_rotset;
		ref_rotset.resize( 20 );  //there are 20 residues in the test protein

		for ( Size i = 1; i <= 20; ++i ) {
			if ( ( i == 8 ) || ( i == 11 ) ) num_rotamers.push_back( 1 );
			else {
				TR << rotsets->nrotamers_for_moltenres( rotsets->resid_2_moltenres( i ) ) << " rotamers at position " << i << std::endl;
				num_rotamers.push_back( rotsets->nrotamers_for_moltenres( rotsets->resid_2_moltenres( i ) ) );
			}
		}

		utility::vector1< pose::Pose > ref_poses;
		core::import_pose::pose_from_file( ref_poses, "core/pack/rotamer_set/test_rotamer_sets.pdb" , core::import_pose::PDB_file);

		for ( Size i = 1; i <= ref_poses.size(); ++i ) {
			for ( Size j = 1; j <= ref_poses[ i ].size(); ++j ) {

				Size pdb_res = ref_poses[ i ].pdb_info()->number( j );
				char pdb_chain = ref_poses[ i ].pdb_info()->chain( j );

				Size this_res = pose.pdb_info()->pdb2pose( pdb_chain, pdb_res );

				ResidueCOP res = ref_poses[ i ].residue( j ).get_self_ptr();
				ref_rotset[ this_res ].push_back( res );
			}
		} //loop over poses

		// AMW: reactivate these tests once the issue with single residue PDB dumping is fixed.
		//now let's go through all the rotamers and compare them from ref to newly generated set
		for ( Size i = 1; i <= 20; ++i ) {

			//TS_ASSERT_EQUALS( ref_rotset[ i ].size(), num_rotamers[ i ] );

			if ( ( i == 8 ) || ( i == 11 ) ) continue;

			Size this_moltenres = rotsets->resid_2_moltenres( i );

			for ( Size j = 1; j <= num_rotamers[ i ]; ++j ) {

				Size this_rotamer_no = rotsets->nrotamer_offset_for_moltenres( this_moltenres ) + j;

				//super stringent: we're making sure that the chi of every residue is correct

				// first make sure that the chi arrays are the same size (ALF notes that sometimes errors occur with database updates --PDR)
				//TS_ASSERT_EQUALS( ref_rotset[ i ][ j ]->nchi(), rotsets->rotamer( this_rotamer_no )->nchi() );
				// cxxtest doesn't seem to provide a macro to test for a condition and then stop the test based on that conditional
				// if the two valuse tested above below are not equal it is best to just return as proceeding any further will cause
				// a vector overrun all tests will stop
				if ( ref_rotset[ i ][ j ]->nchi() != rotsets->rotamer( this_rotamer_no )->nchi() ) {
					return;
				}

				for ( Size k = 1; k <= ref_rotset[ i ][ j ]->nchi(); ++k ) {
					//have to translate chis to periodicity to avoid 179->-180 singularity of dihedrals
					//core::Real ref_chi = basic::periodic_range( rad_per_deg * ref_rotset[ i ][ j ]->chi( k ), pi );
					//core::Real new_chi = basic::periodic_range( rad_per_deg * rotsets->rotamer( this_rotamer_no )->chi( k ), pi );
					//TS_ASSERT_DELTA( ref_chi , new_chi , 0.01 );
				}
			}
		}
		TR << "Done RotamerSets test... " << std::endl;
	}

	void test_rotamer_set_residue_type_offsets() {
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::pack::task;
		using namespace core::pack::rotamer_set;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.8 );
		sfxn.set_weight( fa_rep, 0.44 );
		sfxn( *trpcage );

		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			if ( ii == 5 ) continue;
			task->nonconst_residue_task( ii ).prevent_repacking();
		}

		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( *trpcage, sfxn, task );

		RotamerSet_ res5rots;
		res5rots.set_resid( 5 );
		res5rots.build_rotamers( *trpcage, sfxn, *task, packer_neighbor_graph, true );

		//std::cout << "res5rots.num_rotamers() " << res5rots.num_rotamers() << std::endl;
		Size count_restypes = 1;
		Size count_rots_per_restype = 1;
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );
	}

	Size
	restype_index_for_rt(
		core::pack::rotamer_set::RotamerSet const & rotset,
		std::string const & query_rtname
	)
	{
		for ( Size ii = 1; ii <= rotset.get_n_residue_types(); ++ii ) {
			Size ii_first_rot = rotset.get_residue_type_begin( ii );
			if ( rotset.rotamer_ref( ii_first_rot ).type().name() == query_rtname ) {
				return ii;
			}
		}
		return 0;
	}

	void test_rotamer_set_add_rotamer_into_existing_group1() {
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::pack::task;
		using namespace core::pack::rotamer_set;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.8 );
		sfxn.set_weight( fa_rep, 0.44 );
		sfxn( *trpcage );

		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		task->or_include_current( true );
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			if ( ii == 5 ) continue;
			task->nonconst_residue_task( ii ).prevent_repacking();
		}

		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( *trpcage, sfxn, task );

		RotamerSet_ res5rots;
		res5rots.set_resid( 5 );
		res5rots.build_rotamers( *trpcage, sfxn, *task, packer_neighbor_graph, true );
		Size orig_gln_curr_ind = res5rots.id_for_current_rotamer();
		TS_ASSERT( orig_gln_curr_ind != 0 );
		conformation::ResidueCOP curr_rot = res5rots.rotamer( orig_gln_curr_ind );

		//std::cout << "res5rots.num_rotamers() " << res5rots.num_rotamers() << std::endl;
		Size count_restypes = 1;
		Size count_rots_per_restype = 1;
		Size orig_nrots = res5rots.num_rotamers();
		Size orig_nrts = res5rots.get_n_residue_types();
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );

		Size first_CYS_rot = res5rots.get_residue_type_begin( restype_index_for_rt( res5rots, "CYS" ) );
		Size first_GLN_rot = res5rots.get_residue_type_begin( restype_index_for_rt( res5rots, "GLN" ) );
		Size first_TRP_rot = res5rots.get_residue_type_begin( restype_index_for_rt( res5rots, "TRP" ) );
		//std::cout << "Cys: " << first_CYS_rot << std::endl;
		//std::cout << "Gln: " << first_GLN_rot << std::endl;
		//std::cout << "Trp: " << first_TRP_rot << std::endl;
		TS_ASSERT( first_GLN_rot != 0 && first_TRP_rot != 0 );

		conformation::ResidueCOP rot_CYS = res5rots.rotamer(first_CYS_rot);
		conformation::ResidueCOP rot_GLN = res5rots.rotamer(first_GLN_rot);
		conformation::ResidueCOP rot_TRP = res5rots.rotamer(first_TRP_rot);
		res5rots.add_rotamer_into_existing_group( *rot_GLN );
		res5rots.add_rotamer_into_existing_group( *rot_GLN );
		res5rots.add_rotamer_into_existing_group( *rot_GLN );
		res5rots.add_rotamer_into_existing_group( *rot_TRP );
		res5rots.add_rotamer_into_existing_group( *rot_TRP );
		res5rots.add_rotamer_into_existing_group( *rot_TRP );
		res5rots.add_rotamer_into_existing_group( *rot_TRP );
		res5rots.add_rotamer_into_existing_group( *rot_CYS );
		res5rots.add_rotamer_into_existing_group( *rot_CYS );

		TS_ASSERT_EQUALS( res5rots.get_n_residue_types(), orig_nrts ); // should not change
		TS_ASSERT_EQUALS( res5rots.num_rotamers(), orig_nrots + 9 );
		TS_ASSERT_EQUALS( res5rots.id_for_current_rotamer(), orig_gln_curr_ind + 2 );
		TS_ASSERT_EQUALS( res5rots.rotamer( res5rots.id_for_current_rotamer() ), curr_rot );

		count_rots_per_restype = 1; count_restypes = 1;
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );

	}

	// Variant #2
	// Insert HIS into a PackerTask that does not have HIS already
	void test_rotamer_set_add_rotamer_into_existing_group2() {
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::pack::task;
		using namespace core::pack::rotamer_set;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.8 );
		sfxn.set_weight( fa_rep, 0.44 );
		sfxn( *trpcage );

		PackerTaskOP task1 = TaskFactory::create_packer_task( *trpcage );
		PackerTaskOP task2 = TaskFactory::create_packer_task( *trpcage );

		task1->or_include_current( true );
		utility::vector1< bool > no_HY( chemical::num_canonical_aas, true );
		no_HY[ chemical::aa_his ] = false;
		no_HY[ chemical::aa_tyr ] = false;
		utility::vector1< bool > only_HY( chemical::num_canonical_aas, false );
		only_HY[ chemical::aa_his ] = true;
		only_HY[ chemical::aa_tyr ] = true;

		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			task1->nonconst_residue_task( ii ).restrict_absent_canonical_aas( no_HY );
			task2->nonconst_residue_task( ii ).restrict_absent_canonical_aas( only_HY );
			if ( ii == 5 ) continue;
			task1->nonconst_residue_task( ii ).prevent_repacking();
			task2->nonconst_residue_task( ii ).prevent_repacking();
		}

		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( *trpcage, sfxn, task1 );

		RotamerSet_ res5rots;
		RotamerSet_ res5rots_his;
		res5rots.set_resid( 5 );
		res5rots_his.set_resid( 5 );
		res5rots.build_rotamers( *trpcage, sfxn, *task1, packer_neighbor_graph, true );
		res5rots_his.build_rotamers( *trpcage, sfxn, *task2, packer_neighbor_graph, true );

		Size orig_gln_curr_ind = res5rots.id_for_current_rotamer();
		TS_ASSERT( orig_gln_curr_ind != 0 );
		conformation::ResidueCOP curr_rot = res5rots.rotamer( orig_gln_curr_ind );

		//std::cout << "res5rots.num_rotamers() " << res5rots.num_rotamers() << std::endl;
		Size count_restypes = 1;
		Size count_rots_per_restype = 1;
		Size orig_nrots = res5rots.num_rotamers();
		Size orig_nrts = res5rots.get_n_residue_types();
		Size orig_nrgs = res5rots.get_n_residue_groups();
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );

		Size first_HIS_rot = res5rots_his.get_residue_type_begin( restype_index_for_rt( res5rots_his, "HIS" ) );
		Size first_TYR_rot = res5rots_his.get_residue_type_begin( restype_index_for_rt( res5rots_his, "TYR" ) );

		conformation::ResidueCOP rot_HIS = res5rots_his.rotamer(first_HIS_rot);
		conformation::ResidueCOP rot_TYR = res5rots_his.rotamer(first_TYR_rot);
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );

		TS_ASSERT_EQUALS( res5rots.get_n_residue_types(), orig_nrts+2 );
		TS_ASSERT_EQUALS( res5rots.get_n_residue_types(), orig_nrgs+2 );
		TS_ASSERT_EQUALS( res5rots.num_rotamers(), orig_nrots + 9 );
		TS_ASSERT_EQUALS( res5rots.id_for_current_rotamer(), orig_gln_curr_ind );
		TS_ASSERT_EQUALS( res5rots.rotamer( res5rots.id_for_current_rotamer() ), curr_rot );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "HIS" ), orig_nrts+1 );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "TYR" ), orig_nrts+2 );

		count_rots_per_restype = 1; count_restypes = 1;
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );


		// Now that we have a RotamerSet with HIS but not HIS_D; let's add some HIS_D to make sure
		// that it is put into the same group as HIS.
		Size first_HIS_D_rot = res5rots_his.get_residue_type_begin( restype_index_for_rt( res5rots_his, "HIS_D" ) );
		conformation::ResidueCOP rot_HIS_D = res5rots_his.rotamer(first_HIS_D_rot);
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );

		TS_ASSERT_EQUALS( res5rots.get_n_residue_types(), orig_nrts+3 );
		TS_ASSERT_EQUALS( res5rots.get_n_residue_groups(), orig_nrgs+2 );
		TS_ASSERT_EQUALS( res5rots.num_rotamers(), orig_nrots + 13 );
		TS_ASSERT_EQUALS( res5rots.id_for_current_rotamer(), orig_gln_curr_ind );
		TS_ASSERT_EQUALS( res5rots.rotamer( res5rots.id_for_current_rotamer() ), curr_rot );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "HIS" ), orig_nrts+1 );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "HIS_D" ), orig_nrts+2 );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "TYR" ), orig_nrts+3 );

		count_rots_per_restype = 1; count_restypes = 1;
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );

	}

	// Final variant: add both HIS and HIS_D at the same time.
	void test_rotamer_set_add_rotamer_into_existing_group3() {
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::pack::task;
		using namespace core::pack::rotamer_set;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.8 );
		sfxn.set_weight( fa_rep, 0.44 );
		sfxn( *trpcage );

		PackerTaskOP task1 = TaskFactory::create_packer_task( *trpcage );
		PackerTaskOP task2 = TaskFactory::create_packer_task( *trpcage );

		task1->or_include_current( true );
		utility::vector1< bool > no_HY( chemical::num_canonical_aas, true );
		no_HY[ chemical::aa_his ] = false;
		no_HY[ chemical::aa_tyr ] = false;
		utility::vector1< bool > only_HY( chemical::num_canonical_aas, false );
		only_HY[ chemical::aa_his ] = true;
		only_HY[ chemical::aa_tyr ] = true;

		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			task1->nonconst_residue_task( ii ).restrict_absent_canonical_aas( no_HY );
			task2->nonconst_residue_task( ii ).restrict_absent_canonical_aas( only_HY );
			if ( ii == 5 ) continue;
			task1->nonconst_residue_task( ii ).prevent_repacking();
			task2->nonconst_residue_task( ii ).prevent_repacking();
		}

		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( *trpcage, sfxn, task1 );

		RotamerSet_ res5rots;
		RotamerSet_ res5rots_his;
		res5rots.set_resid( 5 );
		res5rots_his.set_resid( 5 );
		res5rots.build_rotamers( *trpcage, sfxn, *task1, packer_neighbor_graph, true );
		res5rots_his.build_rotamers( *trpcage, sfxn, *task2, packer_neighbor_graph, true );

		Size orig_gln_curr_ind = res5rots.id_for_current_rotamer();
		TS_ASSERT( orig_gln_curr_ind != 0 );
		conformation::ResidueCOP curr_rot = res5rots.rotamer( orig_gln_curr_ind );

		//std::cout << "res5rots.num_rotamers() " << res5rots.num_rotamers() << std::endl;
		Size count_restypes = 1;
		Size count_rots_per_restype = 1;
		Size orig_nrots = res5rots.num_rotamers();
		Size orig_nrts = res5rots.get_n_residue_types();
		Size orig_nrgs = res5rots.get_n_residue_groups();
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );

		Size first_HIS_rot = res5rots_his.get_residue_type_begin( restype_index_for_rt( res5rots_his, "HIS" ) );
		Size first_TYR_rot = res5rots_his.get_residue_type_begin( restype_index_for_rt( res5rots_his, "TYR" ) );
		Size first_HIS_D_rot = res5rots_his.get_residue_type_begin( restype_index_for_rt( res5rots_his, "HIS_D" ) );

		conformation::ResidueCOP rot_HIS = res5rots_his.rotamer(first_HIS_rot);
		conformation::ResidueCOP rot_TYR = res5rots_his.rotamer(first_TYR_rot);
		conformation::ResidueCOP rot_HIS_D = res5rots_his.rotamer(first_HIS_D_rot);
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_HIS );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );
		res5rots.add_rotamer_into_existing_group( *rot_TYR );
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );
		res5rots.add_rotamer_into_existing_group( *rot_HIS_D );

		TS_ASSERT_EQUALS( res5rots.get_n_residue_types(), orig_nrts+3 );
		TS_ASSERT_EQUALS( res5rots.get_n_residue_groups(), orig_nrgs+2 );
		TS_ASSERT_EQUALS( res5rots.num_rotamers(), orig_nrots + 13 );
		TS_ASSERT_EQUALS( res5rots.id_for_current_rotamer(), orig_gln_curr_ind );
		TS_ASSERT_EQUALS( res5rots.rotamer( res5rots.id_for_current_rotamer() ), curr_rot );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "HIS" ), orig_nrts+1 );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "HIS_D" ), orig_nrts+2 );
		TS_ASSERT_EQUALS( restype_index_for_rt( res5rots, "TYR" ), orig_nrts+3 );

		count_rots_per_restype = 1; count_restypes = 1;
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );
	}

	// Variant 4:
	// Add more GLN rotamers to a residue that already is GLN
	void test_rotamer_set_add_rotamer_into_existing_group4() {
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::pack::task;
		using namespace core::pack::rotamer_set;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.8 );
		sfxn.set_weight( fa_rep, 0.44 );
		sfxn( *trpcage );

		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		task->or_include_current( true );
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			if ( ii == 5 ) task->nonconst_residue_task( ii ).restrict_to_repacking();
			else           task->nonconst_residue_task( ii ).prevent_repacking();
		}

		utility::graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( *trpcage, sfxn, task );

		RotamerSet_ res5rots;
		res5rots.set_resid( 5 );
		res5rots.build_rotamers( *trpcage, sfxn, *task, packer_neighbor_graph, true );
		Size orig_gln_curr_ind = res5rots.id_for_current_rotamer();
		TS_ASSERT( orig_gln_curr_ind != 0 );
		conformation::ResidueCOP curr_rot = res5rots.rotamer( orig_gln_curr_ind );

		//std::cout << "res5rots.num_rotamers() " << res5rots.num_rotamers() << std::endl;
		Size count_restypes = 1;
		Size count_rots_per_restype = 1;
		Size orig_nrots = res5rots.num_rotamers();
		Size orig_nrts = res5rots.get_n_residue_types();
		for ( Size ii = 2; ii <= res5rots.num_rotamers(); ++ii ) {
			if ( & res5rots.rotamer(ii)->type() != & res5rots.rotamer(ii-1)->type() ) {
				//std::cout << "nrots for " << res5rots.rotamer(ii-1)->type().name() << ": " << count_rots_per_restype << std::endl;
				TS_ASSERT_EQUALS( count_rots_per_restype, res5rots.get_n_rotamers_for_residue_type( count_restypes ));
				count_rots_per_restype = 1;
				++count_restypes;
				TS_ASSERT_EQUALS( res5rots.get_residue_type_begin( count_restypes ), ii );
			} else {
				++count_rots_per_restype;
			}
			TS_ASSERT_EQUALS( res5rots.get_residue_type_index_for_rotamer( ii ), count_restypes );
		}
		TS_ASSERT_EQUALS( count_restypes, res5rots.get_n_residue_types() );

		Size first_GLN_rot = res5rots.get_residue_type_begin( restype_index_for_rt( res5rots, "GLN" ) );
		//std::cout << "Gln: " << first_GLN_rot << std::endl;
		TS_ASSERT( first_GLN_rot != 0  );

		conformation::ResidueCOP rot_GLN = res5rots.rotamer(first_GLN_rot);
		res5rots.add_rotamer_into_existing_group( *rot_GLN );
		res5rots.add_rotamer_into_existing_group( *rot_GLN );
		res5rots.add_rotamer_into_existing_group( *rot_GLN );

		TS_ASSERT_EQUALS( res5rots.get_n_residue_types(), orig_nrts ); // should not change
		TS_ASSERT_EQUALS( res5rots.num_rotamers(), orig_nrots + 3 );
		TS_ASSERT_EQUALS( res5rots.id_for_current_rotamer(), orig_gln_curr_ind );
		TS_ASSERT_EQUALS( res5rots.rotamer( res5rots.id_for_current_rotamer() ), curr_rot );

	}




};
