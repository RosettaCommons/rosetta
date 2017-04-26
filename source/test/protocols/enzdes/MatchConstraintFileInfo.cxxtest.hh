// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/enzdes/MatchConstraintFileInfo.cxxtest.hh
/// @brief  test suite for constraints between protein and ligand
/// @author Florian Richter

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
//#include <core/conformation/Residue.hh>

//#include <core/io/pdb/pdb_writer.hh>

//#include <core/kinematics/MoveMap.hh>

//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/MinimizerOptions.hh>

//#include <core/pose/Pose.hh>

//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/constraints/ConstraintSet.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/chemical/ChemicalManager.hh> //need for additional residue
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/options/option.hh> //needed to set option
//#include <core/scoring/constraints/AngleConstraint.hh>
//#include <core/scoring/constraints/DihedralConstraint.hh>
//#include <core/scoring/func/Func.hh>
//#include <core/scoring/func/HarmonicFunc.hh>
//#include <core/scoring/constraints/BoundConstraint.hh> //need function in this file
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>


//minimization stuff
//#include <core/kinematics/MoveMap.hh>
//#include <protocols/simple_moves/MinMover.hh>
//#include <protocols/moves/Mover.hh>

#include <core/types.hh>
//#include <math.h>  //need for sqrt taking

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.enzdes.MatchConstraintFileInfo.cxxtest");

using namespace core;


class MatchConstraintFileInfoTest : public CxxTest::TestSuite
{

public:
	MatchConstraintFileInfoTest() {};
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP enz_io;


	// Shared initialization goes here.
	void setUp() {
		using namespace core::chemical;

		core_init_with_additional_options("-run:preserve_header -extra_res_fa protocols/enzdes/D2MX.params");
		/*
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		utility::vector1< std::string > params_files;
		*/
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		/*
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if ( !residue_set.has_name("D2MX") ) params_files.push_back("protocols/enzdes/D2MX.params");
		residue_set.read_files_for_custom_residue_types(params_files);
		basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);
		*/

		enz_io = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO(const_residue_set) );


	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_interface_constraints()
	{
		using namespace protocols::toolbox::match_enzdes_util;

		core::chemical::ResidueTypeSetCOP const_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		//typedef core::id::AtomID AtomID;

		//now let's use the enzdes machinery to read in a cstfile and generate
		//the constraint set, results should be identical to manually created constraints
		enz_io->read_enzyme_cstfile("protocols/enzdes/mcfi_test.cst");

		//enz_io->clear_pdb_specific_data();
		//enz_io->process_pdb_header(compare_pose, catalytic_res);
		//enz_io->check_data_consistency(compare_pose);

		//1. asserting stuff for residue 1
		MatchConstraintFileInfoListCOP mcfil1 = enz_io->mcfi_list( 1 );

		utility::vector1< core::chemical::ResidueTypeCOP > const & up_res1 = mcfil1->upstream_restypes();

		TS_ASSERT_EQUALS( up_res1.size(), 2 );

		// yab 20091021: The following test originally tested explicitly for HIS
		// at up_res1[1] and HIS_D at up_res1[2].  Unfortunately this is a bit
		// problematic as a unit test because the restypes are not loaded into
		// the upstream_restypes container in a consistent order.  What happens
		// in the code is that the restypes are first gathered in a set< CAP >.
		// The set is then iterated over and the restypes loaded in the final
		// vector.  This means the order of the restypes in the vector is
		// effectively dependent on memory addresses, which could be allocated
		// any which way.  Changing this so that we require either HIS - HIS_D
		// or HIS_D - HIS.
		bool const found_his = ( up_res1[ 1 ]->name() == const_residue_set->name_map( "HIS" ).name() ) ||
			( up_res1[ 2 ]->name() == const_residue_set->name_map( "HIS" ).name() );

		bool const found_his_d = ( up_res1[ 1 ]->name() == const_residue_set->name_map( "HIS_D" ).name() ) ||
			( up_res1[ 2 ]->name() == const_residue_set->name_map( "HIS_D" ).name() );

		TS_ASSERT( found_his && found_his_d );

		TS_ASSERT_EQUALS( mcfil1->mcfi(1)->is_covalent(), true );

		//done w stuff for residue 1

		//2. asserting stuff for residue 2
		MatchConstraintFileInfoListCOP mcfil2 = enz_io->mcfi_list( 2 );

		//TO DO


		utility::vector1< std::string > const & algo_strings = mcfil2->mcfi(1)->algorithm_inputs().find( "test" )->second;

		TS_ASSERT_EQUALS( algo_strings.size(), 3 );
		TS_ASSERT_EQUALS( algo_strings[2], "  not so fat after all" );

		utility::vector1< core::Size > const & at_ind_upres2_ser
			= mcfil2->mcfi(1)->enz_cst_template_res( 2 )->atom_inds_for_restype( 2, const_residue_set->name_mapOP("SER") );

		TS_ASSERT_EQUALS( at_ind_upres2_ser.size(), 1 );
		TS_ASSERT_EQUALS( at_ind_upres2_ser[1], const_residue_set->name_map("SER").atom_index("CB")  );

		//done w stuff for residue 2

		//3. asserting stuff for residue 3
		MatchConstraintFileInfoListCOP mcfil3 = enz_io->mcfi_list( 3 );

		TS_ASSERT( mcfil3->mcfi(1)->create_exgs().size() == 0 );

		utility::vector1< SingleConstraint > constraints( mcfil3->mcfi(2)->constraints() );
		TS_ASSERT_EQUALS( constraints.size(), 1 );
		SingleConstraint const & constraint = constraints[ 1 ];

		TS_ASSERT_EQUALS( constraint.tor_U1D3->num_steps(), 4 );
		TS_ASSERT_EQUALS( mcfil3->mcfi(2)->template_atom_inds( 2, 1, const_residue_set->name_map("THR") )[1], const_residue_set->name_map("THR").atom_index("OG1") );
		//TO DO

		//done w stuff for residue 3


		//4. asserting stuff for residue 4
		MatchConstraintFileInfoListCOP mcfil4 = enz_io->mcfi_list( 4 );

		utility::vector1< core::chemical::ResidueTypeCOP > const & up_res4 = mcfil4->upstream_restypes();
		TS_ASSERT_EQUALS( up_res4.size(), 3 );


		utility::vector1< core::Size > const & down_res4_at1 = mcfil4->mcfi(1)->template_atom_inds( 1, 1, const_residue_set->name_map("D2MX") );

		TS_ASSERT_EQUALS( down_res4_at1.size(), 1 );

		TS_ASSERT_EQUALS( down_res4_at1[1], const_residue_set->name_map("D2MX").atom_index("X1") );

		utility::vector1< core::Size > const & at_ind_upres4_phe = mcfil4->mcfi(1)->enz_cst_template_res( 2 )->atom_inds_for_restype( 1, const_residue_set->name_mapOP("PHE") );

		TS_ASSERT_EQUALS( at_ind_upres4_phe.size(), 6 );

		//TS_ASSERT_EQUALS

		//done w stuff for residue 4

		//TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::angle_constraint ] , test_pose.energies().total_energies()[ scoring::angle_constraint ], 1e-5 );
		//TS_ASSERT_DELTA( compare_pose.energies().total_energies()[ scoring::dihedral_constraint ], test_pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-5 );


	}


};
