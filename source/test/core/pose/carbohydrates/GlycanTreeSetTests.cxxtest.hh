// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pose/carbohydrates/GlycanTreeSetTests.cxxtest.hh
/// @brief  Tests for GlycanTreeSet.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Project Headers
#include <core/conformation/carbohydrates/GlycanTreeSetObserver.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>
#include <core/conformation/carbohydrates/GlycanNode.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/carbohydrates/util.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/carbohydrates/util.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("GlycanTreeSetTests");

//using namespace core::pose::carbohydrates;
using namespace core::conformation::carbohydrates;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace core::conformation;
using namespace core::chemical;

class GlycanTreeSetTests : public CxxTest::TestSuite {

public:
	void setUp(){
		core_init_with_additional_options("-include_sugars");
		pose_from_file(pose_, "core/chemical/carbohydrates/gp120_2glycans_man5.pdb", PDB_file);
		pose_from_file( glycan_free_pose_, "core/pose/onechain.pdb", PDB_file);

	}

	void run_basic_tests(core::pose::Pose const & pose, GlycanTreeSetCOP tree_set, utility::vector1< core::Size > const & start_points){

		TR << "Running basic tests" << std::endl;
		//Check that residues downstream are correctly in each tree, and assigned correct values.
		for ( core::Size start_point : start_points ) {
			utility::vector1< core::Size > downstream_residues_1 = get_carbohydrate_residues_of_branch(pose.conformation(), start_point);

			for ( core::Size resnum : downstream_residues_1 ) {
				TS_ASSERT( tree_set->get_tree( start_point )->has_glycan_residue(resnum));
				TS_ASSERT_EQUALS( tree_set->get_tree_containing_residue(resnum)->get_root(), tree_set->get_tree( start_point )->get_root());

				//Assert that each of these are the correct nodes
				TS_ASSERT_EQUALS(tree_set->get_node( resnum )->get_resnum(), resnum );

				//Assert hat the parents are correctly assigned
				TS_ASSERT_EQUALS(tree_set->get_node( resnum )->get_parent(), find_seqpos_of_saccharides_parent_residue(pose.residue( resnum )));

				//Assert that the logic of res to root is correct.
				TS_ASSERT_EQUALS(tree_set->get_tree_root_of_glycan_residue( resnum ), tree_set->get_tree( start_point )->get_root());

				//Assert that we are correctly assigning the exocyclic nature of each residue.
				TS_ASSERT_EQUALS(tree_set->has_exocyclic_glycosidic_linkage( resnum ), has_exocyclic_glycosidic_linkage( pose.conformation(), resnum ) );

				//Assert correct distance to root.
				TS_ASSERT_EQUALS(tree_set->get_distance_to_start( resnum ), get_distance_to_start( pose.conformation(), resnum ) );

				//Assert correct mainchain child.
				TS_ASSERT_EQUALS(tree_set->get_node( resnum )->get_mainchain_child(), find_seqpos_of_saccharides_mainchain_child( pose.residue( resnum )));

				//Make sure linkage position is correct.
				TS_ASSERT_EQUALS(tree_set->get_node( resnum )->get_linkage_position(), get_linkage_position_of_saccharide_residue( pose.residue( resnum ), pose.residue( find_seqpos_of_saccharides_parent_residue(pose.residue( resnum )) ) ) );

			}
		}

	}
	void test_non_glycan_pose() {
		TR << "Testing non-glycan pose" << std::endl;
		TS_ASSERT_EQUALS( glycan_free_pose_.glycan_tree_set(), nullptr);
		TR << "Complete" <<std::endl;

	}
	void test_residue_deletions_and_additions() {
		TR << "Testing residue addition and deletion" << std::endl;
		pose_.delete_polymer_residue(5);


		TR << "residue deleted. " << std::endl;

		GlycanTreeSetCOP tree_set = pose_.glycan_tree_set();
		utility::vector1< GlycanTreeCOP > trees( tree_set->get_all_trees() );


		TS_ASSERT_EQUALS( tree_set->n_trees(), 2);
		TS_ASSERT_EQUALS( tree_set->has_tree( 592 ), false);
		TS_ASSERT_EQUALS( tree_set->has_tree( 585 ), false);

		run_basic_tests( pose_, tree_set, tree_set->get_start_points());

		//Add a residue before the other one we deleted.
		ResidueTypeSetCOP residue_set
			( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		core::conformation::ResidueOP ala_rsd( ResidueFactory::create_residue( residue_set->name_map( "ALA" ) ) );
		pose_.prepend_polymer_residue_before_seqpos(*ala_rsd, 5, true);



		tree_set = pose_.glycan_tree_set();
		trees =  tree_set->get_all_trees() ;


		TS_ASSERT_EQUALS( tree_set->n_trees(), 2);
		TS_ASSERT_EQUALS( tree_set->has_tree( 592 ), true);
		TS_ASSERT_EQUALS( tree_set->has_tree( 585 ), true);

		run_basic_tests( pose_, tree_set, tree_set->get_start_points());

		//Add a pose at the the end of the pose.
		pose_.append_pose_by_jump( glycan_free_pose_, 584);
		tree_set = pose_.glycan_tree_set();

		trees = tree_set->get_all_trees();


		TS_ASSERT_EQUALS( tree_set->n_trees(), 2);
		TS_ASSERT_EQUALS( tree_set->has_tree( 592 ), true);
		TS_ASSERT_EQUALS( tree_set->has_tree( 585 ), true);

		run_basic_tests( pose_, tree_set, tree_set->get_start_points());

	}
	void assert_tree_set_observation_on_end_deletion( core::pose::Pose & pose){
		pose.delete_residue_slow(591);

		GlycanTreeSetCOP tree_set = pose.glycan_tree_set();
		utility::vector1< GlycanTreeCOP > const trees =  tree_set->get_all_trees() ;

		TS_ASSERT_EQUALS(tree_set->n_trees(), 2 );
		TS_ASSERT_EQUALS(tree_set->get_tree(585)->size(), 6);
		TR << "Complete" << std::endl;

	}
	void test_glycan_end_residue_deletion() {
		TR << "Testing end residue deletion" <<std::endl;
		//Test deletion of end of glycan residue.
		assert_tree_set_observation_on_end_deletion( pose_ );
	}

	/* Mid residue deletion is not actually supported yet.
	void
	test_glycan_mid_residue_deletion() {
	TR << "Testing mid residue deletion" << std::endl;
	//Test deletion of mid-glycan.  This will split the glycan trees into multiple.
	pose_.delete_residue_slow( 593 );

	GlycanTreeSetCOP tree_set = pose_.glycan_tree_set();
	utility::vector1< GlycanTreeCOP > const trees =  tree_set->get_all_trees() ;

	TR << "NTrees" << tree_set->size() << std::endl;
	TR << utility::to_string( tree_set->get_start_points()) << std::endl;
	TS_ASSERT_EQUALS( tree_set->n_trees(), 3 );
	TS_ASSERT_EQUALS( tree_set->get_tree( 585 )->size(), 7 );
	TS_ASSERT_EQUALS( tree_set->get_tree( 592 )->size(), 1 );
	TS_ASSERT_EQUALS( tree_set->get_tree( 593 )->size(), 5 );
	TR << "Complete" << std::endl;


	}
	void
	test_glycan_bp_residue_deletion() {
	TR << "Testing branch-point deletion" << std::endl;

	pose_.delete_residue_slow( 594);

	GlycanTreeSetCOP tree_set = pose_.glycan_tree_set();
	utility::vector1< GlycanTreeCOP > const trees = tree_set->get_all_trees();

	TR << "NTrees" << tree_set->size() << std::endl;

	TS_ASSERT_EQUALS( tree_set->size(), 4 );
	TS_ASSERT_EQUALS( tree_set->get_tree( 585 )->size(), 7 );

	TS_ASSERT_EQUALS( tree_set->get_tree( 592 )->size(), 2 );
	TS_ASSERT_EQUALS( tree_set->get_tree( 594 )->size(), 1 ); //Old 595
	TS_ASSERT_EQUALS( tree_set->get_tree( 595 )->size(), 3 ); //Old 596
	TR << "Complete" <<std::endl;

	}
	*/
	void test_glycosylations() {

		TR << "Testing glycosylating pose" << std::endl;
		std::string man5 = "a-D-Manp-(1->3)-[a-D-Manp-(1->3)-[a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-";
		core::pose::carbohydrates::glycosylate_pose( glycan_free_pose_, glycan_free_pose_.pdb_info()->pdb2pose('A', 21), man5);

		//TS_ASSERT( glycan_free_pose_.glycan_tree_set() != nullptr );

		GlycanTreeSetCOP tree_set = glycan_free_pose_.glycan_tree_set();
		utility::vector1< GlycanTreeCOP > trees( tree_set->get_all_trees() );

		TS_ASSERT_EQUALS( tree_set->n_trees(), 1);
		TS_ASSERT_EQUALS( trees[1]->get_size(), 7);

		utility::vector1< core::Size > start_points = tree_set->get_start_points();
		run_basic_tests( glycan_free_pose_, tree_set, start_points );



		core::Size asn_residue = pose_.pdb_info()->pdb2pose('A', 109);

		core::pose::carbohydrates::glycosylate_pose( pose_, asn_residue, man5, false);
		pose_.dump_pdb("gp120_3_glycans.pdb");


		tree_set = pose_.glycan_tree_set();
		trees = tree_set->get_all_trees();
		TS_ASSERT_EQUALS( tree_set->n_trees(), 3);


		TS_ASSERT_EQUALS( trees[1]->get_size(), 7 );
		TS_ASSERT_EQUALS( trees[2]->get_size(), 7 );
		TS_ASSERT_EQUALS( trees[3]->get_size(), 7 );

		//Check on accessing the correct tree.
		TS_ASSERT(tree_set->has_tree(592));
		TS_ASSERT(tree_set->has_tree(585));

		//Check correct root and identification of the trees.
		TS_ASSERT_EQUALS(tree_set->get_tree(592)->get_root(), find_seqpos_of_saccharides_parent_residue(pose_.residue(592)) );
		TS_ASSERT_EQUALS(tree_set->get_tree(585)->get_root(), find_seqpos_of_saccharides_parent_residue(pose_.residue(585)) );

		TS_ASSERT_EQUALS(tree_set->get_tree(592)->get_start(), 592);
		TS_ASSERT_EQUALS(tree_set->get_tree(585)->get_start(), 585);

		run_basic_tests( pose_, tree_set, tree_set->get_start_points() );

		//core::Size start_point = trees[
	}

	void tearDown(){

	}
	void assert_basic_glycan_tree_set_intact(GlycanTreeSetCOP tree_set){

		//Test number of glycan trees
		utility::vector1< GlycanTreeCOP > const trees = tree_set->get_all_trees();
		TS_ASSERT_EQUALS( trees.size(), 2);
		TS_ASSERT_EQUALS( tree_set->n_trees(), 2);

		//Get length of glycan tree - both trees are man5 of 7 residues.
		TS_ASSERT_EQUALS( trees[1]->get_size(), 7);
		TS_ASSERT_EQUALS( trees[2]->get_size(), 7);

		//Check on accessing the correct tree.
		TS_ASSERT(tree_set->has_tree(592));
		TS_ASSERT(tree_set->has_tree(585));

		//Check correct root and identification of the trees.
		TS_ASSERT_EQUALS(tree_set->get_tree(592)->get_root(), find_seqpos_of_saccharides_parent_residue(pose_.residue(592)) );
		TS_ASSERT_EQUALS(tree_set->get_tree(585)->get_root(), find_seqpos_of_saccharides_parent_residue(pose_.residue(585)) );

		TS_ASSERT_EQUALS(tree_set->get_tree(592)->get_start(), 592);
		TS_ASSERT_EQUALS(tree_set->get_tree(585)->get_start(), 585);

		utility::vector1< core::Size > start_points;
		start_points.push_back( 592 );
		start_points.push_back( 585 );

		run_basic_tests( pose_, tree_set, start_points);
		TR << "Complete" << std::endl;


	}
	void test_glycan_tree_set() {
		TR << "Testing GlycanTreeSet" << std::endl;
		TS_ASSERT( pose_.glycan_tree_set() != nullptr );
		GlycanTreeSetCOP tree_set = pose_.glycan_tree_set();

		assert_basic_glycan_tree_set_intact( tree_set );
	}

	//For now, we are commenting the serialization test out until Andrew gets back.
	// We need this into master ASAP.
	void test_serialization_and_clone() {
		using namespace core::conformation;

		ConformationOP conf = pose_.conformation().clone();
		assert_basic_glycan_tree_set_intact( conf->glycan_tree_set() );

#ifdef SERIALIZATION


		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( conf );
		}

		std::istringstream iss( oss.str() );
		ConformationOP conf2;
		{
			cereal::BinaryInputArchive arc( iss );
			arc( conf2 );
		}

		TR.Debug << "Serialized Form: \n";
		TR.Debug << oss.str() << std::endl;
		assert_basic_glycan_tree_set_intact( conf->glycan_tree_set() );
#endif

	}
	void test_serialization_observation() {
		TS_ASSERT(true);

#ifdef SERIALIZATION

		//Serialze a pose, which SHOULD serialize the underlying Conformation.
		// Test to make sure that the GlycanTreeSet is still observing that conformation.
		core::pose::PoseOP pose = pose_.clone();
		assert_basic_glycan_tree_set_intact( pose->conformation().glycan_tree_set() );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( pose );
		}

		std::istringstream iss( oss.str() );
		core::pose::PoseOP pose2;
		{
			cereal::BinaryInputArchive arc( iss );
			arc( pose2 );
		}

		TR.Debug << "Serialized Form: \n";
		TR.Debug << oss.str() << std::endl;

		assert_tree_set_observation_on_end_deletion( *pose2);

#endif
	}

	// Confirm that Pose::glycan_tree_sequence() returns the correct sequence.
	void test_glycan_tree_sequence()
	{
		using namespace std;
		using namespace core::pose;

		// This is the standard N-linked core of a hybrid, bisected, fucosylated glycan.
		// It adequately tests residues with multiple branches as well as branches off of branches.
		// It also includes modified sugars.
		string const input_sequence(
			"alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)]-"
			"[beta-D-GlcpNAc-(1->4)]-beta-D-Manp-(1->4)-beta-D-GlcpNAc-(1->4)-[alpha-L-Fucp-(1->6)]-beta-D-GlcpNAc" );

		Pose pose;

		make_pose_from_saccharide_sequence( pose, input_sequence );
		TS_ASSERT_EQUALS( pose.glycan_tree_sequence( 1 ), input_sequence );
		TS_ASSERT_EQUALS( pose.glycan_tree_sequence( 2 ), input_sequence );

		make_pose_from_sequence( pose, "ANASA", "fa_standard" );
		core::pose::carbohydrates::glycosylate_pose( pose, 2, input_sequence );
		TS_ASSERT_EQUALS( pose.glycan_tree_sequence( 1 ), "" );
		TS_ASSERT_EQUALS( pose.glycan_tree_sequence( 2 ), "" );
		TS_ASSERT_EQUALS( pose.glycan_tree_sequence( 6 ), input_sequence + "-" );
	}

private:
	core::pose::Pose pose_;
	core::pose::Pose glycan_free_pose_;

};

