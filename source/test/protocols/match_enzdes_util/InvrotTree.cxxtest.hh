// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match_enzdes_util/InvrotTree.cxxtest.hh
/// @brief  test suite for InvrotTreeFunctionality
/// @details this primarily tests the constraints, the inverse rotamer
/// building functionality is covered by cstfile_to_theozyme_pdb
/// integration test
/// @author Florian Richter

// Test headers
#include <cxxtest/TestSuite.h>
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>

#include <test/protocols/init_util.hh>

#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh> //need for additional residue
#include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh> //needed to set option
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.match_enzdes_util.InvrotTree.cxxtest");

using namespace core;

/// @details several things are getting tested
/// 1. simple invrot tree with no protein-protein
/// interactions for 3-residue tim reaction, i.e.
/// the constraint should consist of three ambiguous
/// constraints.
/// 2. more complicated: chorismate mutase or esterase,
/// with some protein-protein interactions. not implemented yet
class InvrotTreeTest : public CxxTest::TestSuite
{

public:
  InvrotTreeTest() {};
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP tim_enz_io;


	// Shared initialization goes here.
	void setUp() {
		protocols_init();
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCAP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if(!residue_set.has_name("1n1")) params_files.push_back("protocols/match_enzdes_util/1n1.params");
		residue_set.read_files(params_files);
		basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);

		tim_enz_io = new protocols::toolbox::match_enzdes_util::EnzConstraintIO(& residue_set);
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  void test_invrot_tree()
  {
	  using namespace core::scoring::constraints;
		using namespace protocols::toolbox::match_enzdes_util;

	  pose::Pose tim_pose;

	  core::import_pose::pose_from_pdb( tim_pose, "protocols/match_enzdes_util/1ney_invtree_unittest.pdb");
	  scoring::ScoreFunctionOP scorefxn = scoring::ScoreFunctionFactory::create_score_function("enzdes");
		scorefxn->set_weight( core::scoring::backbone_stub_constraint, 1.0 );
		(*scorefxn)(tim_pose);
		TR << "At first, pose has coordinate constraint score " << tim_pose.energies().total_energies()[ core::scoring::coordinate_constraint ] << " and backbone stub constraint score " << tim_pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

		// ======= testing simple TIM theozyme ==================
		TR << "starting tim theozyme tree set of tests... " << std::endl;
		//let's build up random attachment pos
		utility::vector1< utility::vector1< Size > > tim_pos(3);
		tim_pos[1].push_back(20);
		//note: 11 in test pose is gly, so that's not going to count
		tim_pos[2].push_back(10); tim_pos[2].push_back(11); tim_pos[2].push_back(12); tim_pos[2].push_back(13);
		tim_pos[3].push_back(4); tim_pos[3].push_back(5);
		AllowedSeqposForGeomCstOP tim_actpos = new AllowedSeqposForGeomCst( tim_pos );


	  tim_enz_io->read_enzyme_cstfile("protocols/match_enzdes_util/mocktim.cst");

		InvrotTreeOP tim_invrot_tree = new TheozymeInvrotTree( tim_enz_io );
		tim_invrot_tree->generate_targets_and_inverse_rotamers();
		utility::vector1< InvrotCollectorCOP > tim_invrots( tim_invrot_tree->collect_all_inverse_rotamers() );
		//tree should contain one unique definition
		TS_ASSERT( tim_invrots.size() == 1 );

		//we have to add the ligand to the pose for constraint generation to work...
		core::conformation::ResidueCOP ligres( *(tim_invrots[1]->invrots()[0].begin()) );
		tim_pose.append_residue_by_jump( *ligres, tim_pose.total_residue() );

		tim_invrot_tree->generate_inverse_rotamer_constraints( tim_pose, tim_actpos );
		TS_ASSERT( tim_invrot_tree->num_target_states() == 1 );

		core::scoring::constraints::ConstraintCOP tim_tree_cst = tim_invrot_tree->get_constraint_for_target_state( 1 );

		//round one of tests: let's examine the constraint a bit
		//should be a multicst containing three ambiguous csts
		TS_ASSERT( dynamic_cast< MultiConstraint const * > ( &(*tim_tree_cst) ) );
		TS_ASSERT( !dynamic_cast< AmbiguousConstraint const * > ( &(*tim_tree_cst) ) );

		MultiConstraint const & tim_tree_multicst( static_cast<MultiConstraint const & > (*tim_tree_cst ) );
		TS_ASSERT( tim_tree_multicst.member_constraints().size() == 3 );
		TS_ASSERT( dynamic_cast< AmbiguousConstraint const * > ( &(*(tim_tree_multicst.member_constraints()[1])) ) );
		TS_ASSERT( dynamic_cast< AmbiguousConstraint const * > ( &(*(tim_tree_multicst.member_constraints()[2])) ) );
		TS_ASSERT( dynamic_cast< AmbiguousConstraint const * > ( &(*(tim_tree_multicst.member_constraints()[3])) ) );

		AmbiguousConstraint const & tim_ambig_cst1( dynamic_cast< AmbiguousConstraint const & > ( *(tim_tree_multicst.member_constraints()[1]) ) );
		AmbiguousConstraint const & tim_ambig_cst2( dynamic_cast< AmbiguousConstraint const & > ( *(tim_tree_multicst.member_constraints()[2]) ) );
		AmbiguousConstraint const & tim_ambig_cst3( dynamic_cast< AmbiguousConstraint const & > ( *(tim_tree_multicst.member_constraints()[3]) ) );

		//check whether the ambig constraints have the expected number of member constraints
		TS_ASSERT( tim_ambig_cst1.member_constraints().size() == tim_invrots[1]->invrots()[1].size() );
		TS_ASSERT( tim_ambig_cst2.member_constraints().size() == (tim_invrots[1]->invrots()[2].size() * 3) );
		TS_ASSERT( tim_ambig_cst3.member_constraints().size() == (tim_invrots[1]->invrots()[3].size() * 2) );

		//now let's add it to the pose test pose and see how it behaves
		tim_pose.add_constraint( tim_tree_cst );
		(*scorefxn)(tim_pose);
		TR << "After adding invrot constraint, pose has coordinate constraint score " << tim_pose.energies().total_energies()[ core::scoring::coordinate_constraint ] << " and backbone stub constraint score " << tim_pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;
		TS_ASSERT_DELTA(1193.3896,  tim_pose.energies().total_energies()[ core::scoring::coordinate_constraint ], 1e-2 );
		TS_ASSERT_DELTA(0.0,  tim_pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ], 1e-2 );

		//now let's set the backbone position in the pose residues to a random
		//one observed in the invrots
		std::list< core::conformation::ResidueCOP >::const_iterator res1_it( tim_invrots[1]->invrots()[1].begin() ), res2_it( tim_invrots[1]->invrots()[2].begin() ), res3_it( tim_invrots[1]->invrots()[3].begin() );

		res1_it++;
		res2_it++; res2_it++;

		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(20).atom_index("N"), 20 ), (*res1_it)->xyz("N") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(20).atom_index("CA"), 20 ), (*res1_it)->xyz("CA") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(20).atom_index("C"), 20 ), (*res1_it)->xyz("C") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(20).atom_index("O"), 20 ), (*res1_it)->xyz("O") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(20).atom_index("CB"), 20 ), (*res1_it)->xyz("CB") );


		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(10).atom_index("N"), 10 ), (*res2_it)->xyz("N") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(10).atom_index("CA"), 10 ), (*res2_it)->xyz("CA") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(10).atom_index("C"), 10 ), (*res2_it)->xyz("C") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(10).atom_index("O"), 10 ), (*res2_it)->xyz("O") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(10).atom_index("CB"), 10 ), (*res2_it)->xyz("CB") );

		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(5).atom_index("N"), 5 ), (*res3_it)->xyz("N") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(5).atom_index("CA"), 5 ), (*res3_it)->xyz("CA") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(5).atom_index("C"), 5 ), (*res3_it)->xyz("C") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(5).atom_index("O"), 5 ), (*res3_it)->xyz("O") );
		tim_pose.set_xyz( core::id::AtomID( tim_pose.residue(5).atom_index("CB"), 5 ), (*res3_it)->xyz("CB") );

		(*scorefxn)(tim_pose);

		TR << "After setting pose res to invrot coords, pose has coordinate constraint score " << tim_pose.energies().total_energies()[ core::scoring::coordinate_constraint ] << " and backbone stub constraint score " << tim_pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;
		TS_ASSERT_DELTA(0.0,  tim_pose.energies().total_energies()[ core::scoring::coordinate_constraint ], 1e-2 );
		TS_ASSERT_DELTA(-60.0,  tim_pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ], 1e-2 );
		//tim_pose.dump_pdb( "tempinvrot_debug.pdb");

		// ======= testing simple TIM theozyme over ==================
	}
};

