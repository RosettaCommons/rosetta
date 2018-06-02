// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/constraint_generator/DihedralConstraintGeneratorTests.cxxtest.hh
/// @brief  Unit tests for the DihedralConstraintGenerator.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <protocols/constraint_generator/DihedralConstraintGenerator.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/residue_selector/CDRResidueSelector.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/types.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <numeric/util.hh>


static basic::Tracer TR("DihedralConstraintGeneratorTests");

using namespace core::pose;
using namespace core::select::residue_selector;
using namespace protocols::constraint_generator;
using namespace core::scoring::constraints;
using namespace core::scoring::func;
using namespace protocols::antibody;
using namespace protocols::antibody::residue_selector;
using namespace core::id;

class DihedralConstraintGeneratorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	core::pose::Pose pose_;
	AntibodyInfoOP ab_info_;

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose_, "protocols/antibody/1bln_AB_aho.pdb", core::import_pose::PDB_file);

		ab_info_ = AntibodyInfoOP(new AntibodyInfo(pose_, AHO_Scheme, North));
	}

	void tearDown(){

	}

	void test_generator(){

		utility::vector1< CDRNameEnum > cdrs;
		cdrs.push_back( l1 ); //Set constraints for L1 residues!
		CDRResidueSelectorOP selector = CDRResidueSelectorOP( new CDRResidueSelector( ab_info_, cdrs ));
		DihedralConstraintGeneratorOP generator = DihedralConstraintGeneratorOP( new DihedralConstraintGenerator());

		generator->set_torsion_type( core::id::phi_dihedral);
		generator->set_residue_selector( selector );

		TS_ASSERT_THROWS_NOTHING( generator->apply( pose_ ) ); //Super basic test just to make sure 'it works'.  Now we do the meat.


		//Test that a single PHI residue is set.  Here, we use 24L since I know it exists in the pose as it is the beginning of the CDR residue.

		core::Size L1_start = ab_info_->get_CDR_start( l1, pose_ );

		utility::vector1< bool > subset(pose_.total_residue(), false);
		subset[L1_start] = true;

		ReturnResidueSubsetSelectorOP subset_selector = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector( subset ));

		generator->set_residue_selector( subset_selector );

		utility::vector1< ConstraintCOP > auto_single_cst = generator->apply( pose_ );
		TS_ASSERT( auto_single_cst.size() == 1); //Assert that this is a constraint for a single residue.


		//Now, we do everything manually and make sure they are the same constraint.


		AtomID C_0 =   AtomID(pose_.residue(L1_start-1).atom_index("C"), L1_start-1);
		AtomID N_1 =   AtomID(pose_.residue(L1_start).atom_index("N"), L1_start );
		AtomID CA_1 = AtomID(pose_.residue(L1_start).atom_index("CA"), L1_start );
		AtomID C_1 =   AtomID(pose_.residue(L1_start).atom_index("C"), L1_start );
		//AtomID N_2 =   AtomID(pose_.residue(L1_start+1).atom_index("N"), L1_start+1 );

		utility::vector1< AtomID > custom_dihedral;
		custom_dihedral.push_back(C_0);
		custom_dihedral.push_back(N_1);
		custom_dihedral.push_back(CA_1);
		custom_dihedral.push_back(C_1);

		generator->set_custom_dihedral( custom_dihedral );

		utility::vector1< ConstraintCOP > manual_single_cst = generator->apply( pose_ );

		ConstraintSet cst_set_auto;
		ConstraintSet cst_set_manual;
		ConstraintSet cst_set_parsed;

		cst_set_auto.add_constraints( auto_single_cst );
		cst_set_manual.add_constraints( manual_single_cst );

		TR << "Auto Constraints" << std::endl;
		cst_set_auto.show_definition(TR, pose_ );
		TR << std::endl << std::endl;

		TR << "Manual Constraints" << std::endl;
		TR << "Manual Constraint Size " << manual_single_cst.size() << std::endl;
		cst_set_manual.show_definition(TR, pose_ );

		TS_ASSERT( manual_single_cst.size() == 1);


		TR << "Auto   " << auto_single_cst[1]->score( pose_ ) << std::endl;
		TR << "Manual " << manual_single_cst[1]->score( pose_ ) << std::endl;


		TS_ASSERT( auto_single_cst[1]->score_type() == manual_single_cst[1]->score_type() );



		TS_ASSERT_DELTA( auto_single_cst[1]->score( pose_ ), manual_single_cst[1]->score( pose_ ), .000001); //Make sure constraints are functionally identical as they should be.


		//Now, we test correct RosettaScripts Parsing.

		std::string atom_names="C,N,CA,C";
		std::string resnums="23L,24L,24L,24L";

		std::string dih = "phi";
		core::id::MainchainTorsionType torsion_type = parse_torsion_type( dih );
		TS_ASSERT( phi_dihedral == torsion_type  );

		utility::vector1< AtomID > parsed_dihedral = parse_custom_torsion( atom_names, resnums, pose_ );

		generator->set_custom_dihedral( parsed_dihedral );
		utility::vector1< ConstraintCOP > parsed_cst = generator->apply( pose_ );

		cst_set_parsed.add_constraints( parsed_cst );
		cst_set_parsed.show_definition( TR, pose_ );
		TS_ASSERT( parsed_cst.size() == 1 );

		TS_ASSERT_DELTA(manual_single_cst[1]->score( pose_ ), parsed_cst[1]->score( pose_ ), .000001);




	}

};
