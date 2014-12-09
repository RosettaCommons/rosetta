// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/DihedralConstraint.cxxtest.hh
/// @brief  test suite for constraints between protein and ligand
/// @author Florian Richter

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/protocols/init_util.hh>

#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh> //need for additional residue
#include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh> //needed to set option

#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh> //need function in this file
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
// AUTO-REMOVED #include <protocols/enzdes/AddorRemoveCsts.hh> //for parser testing
#include <protocols/moves/MoverFactory.hh> //for parser testing

#include <basic/datacache/DataMap.hh> //for parser test
#include <protocols/filters/Filter.hh> //for parser test


//minimization stuff
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/Mover.hh>

// AUTO-REMOVED #include <math.h>  //need for sqrt taking

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <numeric/constants.hh>

//utility headers
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.LigInterfaceConstraints.cxxtest");

using namespace core;


class LigInterfaceConstraintsTest : public CxxTest::TestSuite
{

public:
  LigInterfaceConstraintsTest() {};
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP enz_io;


	// Shared initialization goes here.
	void setUp() {
		protocols_init();
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if(!residue_set.has_name("D2N")) params_files.push_back("protocols/enzdes/D2N.params");
		residue_set.read_files(params_files);
		basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);

		enz_io = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO(residue_set.get_self_weak_ptr()) );
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  void test_interface_constraints()
  {
	  using namespace core::scoring::constraints;
	  typedef core::id::AtomID AtomID;

	  pose::Pose test_pose, compare_pose, parser_compare_pose;
	  core::Real const rad_per_deg = numeric::constants::f::degrees_to_radians;

	  core::import_pose::pose_from_pdb( test_pose, "protocols/enzdes/ligtest_it.pdb");
	  scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
	  scorefxn->reset();
	  //scorefxn->set_weight( scoring::fa_atr, 1.0);
	  scorefxn->set_weight( scoring::atom_pair_constraint, 1.0);
	  scorefxn->set_weight( scoring::angle_constraint, 1.0);
	  scorefxn->set_weight( scoring::dihedral_constraint, 1.0);

	  ConstraintSetOP test_cst_set( new ConstraintSet() );


	  //couple of constraints
	  core::scoring::func::FuncOP his_lig_dist( new core::scoring::constraints::BoundFunc((2.00-0.30),(2.00+0.30), sqrt(1.0/180.0), "dis") );
	  test_cst_set->add_constraint( ConstraintCOP( new AtomPairConstraint( AtomID(test_pose.residue_type(107).atom_index("C6"),107), AtomID(test_pose.residue_type(45).atom_index("ND1"),45), his_lig_dist) ) );
	  TR << "Adding constraint between res " << test_pose.residue_type(45).name3() << "45 and " << test_pose.residue_type(107).name3() <<"107." << std::endl;

	  core::scoring::func::FuncOP his_lig_ang( new core::scoring::constraints::BoundFunc( (105.10*rad_per_deg - 6.00*rad_per_deg), (105.10*rad_per_deg + 6.00*rad_per_deg), sqrt(1.0/100.0),"angA") );
	  test_cst_set->add_constraint( ConstraintCOP( new AngleConstraint( AtomID(test_pose.residue_type(107).atom_index("O4"),107), AtomID(test_pose.residue_type(107).atom_index("C6"),107), AtomID(test_pose.residue_type(45).atom_index("ND1"),45), his_lig_ang) ) );


	  core::scoring::func::FuncOP ser_lig_ang( new core::scoring::constraints::BoundFunc( (109.00*rad_per_deg - 15.00*rad_per_deg), (109.00*rad_per_deg + 15.00*rad_per_deg), sqrt(1.0/20.0),"angB") );
	  test_cst_set->add_constraint( ConstraintCOP( new AngleConstraint( AtomID(test_pose.residue_type(107).atom_index("O4"),107), AtomID(test_pose.residue_type(68).atom_index("OG"),68), AtomID(test_pose.residue_type(68).atom_index("CB"),68), ser_lig_ang) ) );

	  core::scoring::func::FuncOP gln_lig_dist( new core::scoring::constraints::BoundFunc( (3.00-0.20),(3.00 + 0.20), sqrt(1.0/20.0), "dis") );
	  test_cst_set->add_constraint( ConstraintCOP( new AtomPairConstraint( AtomID(test_pose.residue_type(107).atom_index("O4"),107), AtomID(test_pose.residue_type(51).atom_index("NE2"),51), gln_lig_dist) ) );


	  core::scoring::func::FuncOP gln_lig_dih( new core::scoring::constraints::PeriodicBoundFunc( (180.00*rad_per_deg - 15.00*rad_per_deg),(180.00*rad_per_deg + 15.00*rad_per_deg), sqrt(1.0/25.0), "dis", 180.00*rad_per_deg) );
	  test_cst_set->add_constraint( ConstraintCOP( new DihedralConstraint( AtomID(test_pose.residue_type(107).atom_index("O4"),107), AtomID(test_pose.residue_type(51).atom_index("NE2"),51), AtomID(test_pose.residue_type(51).atom_index("CD"),51), AtomID(test_pose.residue_type(51).atom_index("CG"),51), gln_lig_dih) ) );
	  // all constraints defined

	  //test_cst_set->show(TR);

	  test_pose.constraint_set( test_cst_set);
	  (*scorefxn)(test_pose);

	  TR << "Atom Pair Constraint Energy should be 59.3895, and is: " << test_pose.energies().total_energies()[ scoring::atom_pair_constraint ] << ", fa_atr should be 0 and is " << test_pose.energies().total_energies()[ scoring::fa_atr ] << std::endl;

	  TR << "angle Constraint Energy should be 1.67883, and is: " << test_pose.energies().total_energies()[ scoring::angle_constraint ] << std::endl;

	  TR << "dihedral Constraint Energy should be 18.963, and is: " << test_pose.energies().total_energies()[ scoring::dihedral_constraint ] << std::endl;



	  TS_ASSERT_DELTA(59.3895, test_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 1e-4 );
	  TS_ASSERT_DELTA(1.67883, test_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-4 );
	  TS_ASSERT_DELTA(18.963, test_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-3 );


	  //now let's use the enzdes machinery to read in a cstfile and generate
	  //the constraint set, results should be identical to manually created constraints
	  enz_io->read_enzyme_cstfile("protocols/enzdes/ligtest_it.cst");
	  core::import_pose::pose_from_pdb( compare_pose, "protocols/enzdes/ligtest_it.pdb");
		core::pose::Pose parser_pose = compare_pose;


		(*scorefxn)(compare_pose);
		(*scorefxn)(parser_pose);
	  //enz_io->clear_pdb_specific_data();
	  //enz_io->process_pdb_header(compare_pose, catalytic_res);
	  //enz_io->check_data_consistency(compare_pose);
	  enz_io->add_constraints_to_pose(compare_pose, scorefxn, false);
	  (*scorefxn)(compare_pose);

	  TR << "comparing scores of explicit constraint set with cstfile processing generated constraint set... " << std::endl;

	  TR << "atom pair sum for cstfile reading is " << compare_pose.energies().total_energies()[ scoring::atom_pair_constraint ] << std::endl;

	  TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::atom_pair_constraint ] , test_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 1e-6 );

	  TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::angle_constraint ] , test_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-6 );
	  TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::dihedral_constraint ], test_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-6 );


		//we will also test if we can get the same constraint score through using the parser
		TR << "Beginning parser testing setup....   ";
		///protocols::moves::MoverFactory mover_factory; //this isn't static yet
		///mover_factory.add_type( new protocols::enzdes::AddOrRemoveMatchCsts );

		utility::io::izstream fin;
		fin.open("protocols/enzdes/parse_cst_test.xml");
		runtime_assert( fin.good() );
		utility::tag::TagCOP tag = utility::tag::Tag::create( fin );
		fin.close();
		utility::vector0< utility::tag::TagCOP > const TO_tags( tag->getTag("MOVERS")->getTags() );

		protocols::moves::Movers_map parser_movers;
		protocols::filters::Filters_map parser_filters;
		basic::datacache::DataMap parser_datamap; // abstract objects, such as scorefunctions, to be used by filter and movers

		for( utility::vector0< utility::tag::TagCOP >::const_iterator tp( TO_tags.begin() ), tp_e( TO_tags.end() ); tp != tp_e; ++tp ) {

			std::string const user_defined_name( (*tp)->getOption<std::string>("name") );
			protocols::moves::MoverOP new_mover(  protocols::moves::MoverFactory::get_instance()->newMover( *tp, parser_datamap, parser_filters, parser_movers, parser_pose ) );

			parser_movers.insert( std::make_pair( user_defined_name, new_mover ) );
		}

		protocols::moves::MoverOP newcstmover = parser_movers.find("cstaddnew")->second;
		runtime_assert( newcstmover != 0 );
		TR << "parser testing setup done, beginning testing....   " << std::endl;

		newcstmover->apply( parser_pose );
		(*scorefxn)(parser_pose);
		//first testing if adding new works
		TR << "atom pair sum for parser based cstfile reading is " << parser_pose.energies().total_energies()[ scoring::atom_pair_constraint ] << std::endl;
		TS_ASSERT_DELTA( test_pose.energies().total_energies()[ scoring::atom_pair_constraint ] , parser_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 1e-6 );
	  TS_ASSERT_DELTA( parser_pose.energies().total_energies()[ scoring::angle_constraint ] , test_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-6 );
	  TS_ASSERT_DELTA( parser_pose.energies().total_energies()[ scoring::dihedral_constraint ], test_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-6 );

		//then testing if removing works
		protocols::moves::MoverOP remcstmover = parser_movers.find("cstremove")->second;
		runtime_assert( remcstmover != 0 );
		remcstmover->apply( parser_pose );
		(*scorefxn)(parser_pose);

		TR << "atom pair sum for parser based cstfile reading after removing of constraints is " << parser_pose.energies().total_energies()[ scoring::atom_pair_constraint ] << std::endl;
		TS_ASSERT_DELTA( parser_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 0, 1e-6 );
	  TS_ASSERT_DELTA( parser_pose.energies().total_energies()[ scoring::angle_constraint ] , 0, 1e-6 );
	  TS_ASSERT_DELTA( parser_pose.energies().total_energies()[ scoring::dihedral_constraint ], 0, 1e-6 )

		//removing testing done, now testing whether adding pregenerated csts works
		protocols::moves::MoverOP pregencstmover = parser_movers.find("cstaddpreg")->second;
		runtime_assert( pregencstmover != 0 );
		pregencstmover->apply( parser_pose );
		(*scorefxn)(parser_pose);

		TR << "atom pair sum for parser based cstfile reading after readding pregenerated constraints is " << parser_pose.energies().total_energies()[ scoring::atom_pair_constraint ] << std::endl;
		TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::atom_pair_constraint ] , parser_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 1e-6 );
	  TS_ASSERT_DELTA( parser_pose.energies().total_energies()[ scoring::angle_constraint ] , compare_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-6 );
	  TS_ASSERT_DELTA( parser_pose.energies().total_energies()[ scoring::dihedral_constraint ], compare_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-6 );
		//parser testing done


	  //ok, constraint scoring seems to work, now let's do a minimization, shall we?
	  //once again we'll compare results of the explicit constraint set to the one
	  //generated by the cstfile reader

	  Size jump_id = test_pose.num_jump(); // assume ligand attached by last jump
	  TS_ASSERT_EQUALS( test_pose.num_jump(), compare_pose.num_jump());

	  core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	  movemap->set_jump(jump_id, true);
	  movemap->set_chi(45, true);
	  movemap->set_chi(51, true);
	  movemap->set_chi(68, true);


	  protocols::simple_moves::MinMoverOP dfpMinTightTol( new protocols::simple_moves::MinMover( movemap, scorefxn, "dfpmin_armijo_nonmonotone_atol", 0.02, true /*use_nblist*/ ) );

	  protocols::simple_moves::MinMoverOP dfpMinTightTol2( new protocols::simple_moves::MinMover( movemap, scorefxn, "dfpmin_armijo_nonmonotone_atol", 0.02, true /*use_nblist*/ ) );


	  TR << "scoring seems to work, doing minimization...  ";
	  dfpMinTightTol->apply(test_pose);
	  (*scorefxn)(test_pose);

	  dfpMinTightTol2->apply(compare_pose);
	  (*scorefxn)(compare_pose);
	  //		scorefxn->accumulate_residue_total_energies( test_pose );

	  TR << "  done." << std::endl;

	  //and now let's hope that all the constraints have been minimized to equal values that are identical to the precalculated ones

	  TS_ASSERT_DELTA(12.3785, test_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 1e-3 );
	  TS_ASSERT_DELTA(0.3640, test_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-3 );
	  TS_ASSERT_DELTA(0.0051, test_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-3 );

		// TS_ASSERT_DELTA(11.397, test_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 1e-3 );
	  // TS_ASSERT_DELTA(0.168, test_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-3 );
	  // TS_ASSERT_DELTA(0.001, test_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-3 );

	  TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::atom_pair_constraint ] , test_pose.energies().total_energies()[ scoring::atom_pair_constraint ], 1e-5 );

	  TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::angle_constraint ] , test_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-5 );
	  TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::dihedral_constraint ], test_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-5 );



  }


};
