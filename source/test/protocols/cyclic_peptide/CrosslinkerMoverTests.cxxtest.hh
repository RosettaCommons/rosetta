// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/CrosslinkerMoverTests.cxxtest.hh
/// @brief  Unit tests for the CrosslinkerMover.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/cyclic_peptide/CrosslinkerMover.hh>
#include <protocols/cyclic_peptide/crosslinker/TMA_Helper.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/simple_moves/MinMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("CrosslinkerMoverTests");


class CrosslinkerMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options("-extra_res_fa crosslinker/trimesic_acid.params");
	}

	void tearDown(){

	}



	void test_tma_dihedral_fxn(){
		core::pose::Pose pose;
		protocols::cyclic_peptide::PeptideStubMover pepstub;
		pepstub.set_reset_mode(true);
		pepstub.add_residue("Append","ORN:SidechainConjugation:NtermProteinFull:CtermProteinFull",1,true,"",0,0,"");
		pepstub.add_residue("Append","TMA",2,false,"CM1",0,1,"NE");
		pepstub.add_residue("Append", "ORN:SidechainConjugation:NtermProteinFull:CtermProteinFull",3,false,"NE",0,2,"CM2");
		pepstub.add_residue("Append", "ORN:SidechainConjugation:NtermProteinFull:CtermProteinFull",4,false,"NE",0,2,"CM3");
		pepstub.apply(pose);

		//pose.dump_pdb("vtemp1.pdb"); //DELETE ME

		for ( core::Size ir=1; ir<=4; ++ir ) {
			for ( core::Size ia=1, iamax=pose.residue_type(ir).nchi(); ia<=iamax; ++ia ) {
				pose.set_chi(ia,ir,180.0);
			}
		}

		pose.update_residue_neighbors();
		//pose.dump_pdb("vtemp2.pdb"); //DELETE ME

		protocols::cyclic_peptide::crosslinker::TMA_Helper helper;

		core::select::residue_selector::ResidueSubset selection(4, true);
		selection[2] = false;
		helper.add_linker_constraints_asymmetric(pose, selection);

		core::scoring::ScoreFunctionOP sfxn(new core::scoring::ScoreFunction);
		sfxn->set_weight(core::scoring::dihedral_constraint, 1.0);

		TR << "\nChi\tE\n";
		for ( core::Real c(0.0); c<360.1; c+=1.0 ) {
			pose.set_chi(1,2,c);
			pose.set_chi(2,2,c);
			pose.set_chi(3,2,c);
			core::Real const curenergy( (*sfxn)(pose) );
			TR << c << "\t" << curenergy << "\n";
		}
		TR.flush();

		core::kinematics::MoveMapOP movemap(new core::kinematics::MoveMap);
		movemap->set_bb(false);
		movemap->set_jump(false);
		movemap->set_chi(false);
		movemap->set_chi(2, true);
		protocols::simple_moves::MinMover minmove(movemap, sfxn, "dfpmin", 0.0000001, false, false, false);
		pose.set_chi(1,2,30.0);
		pose.set_chi(2,2,30.0);
		pose.set_chi(3,2,30.0);
		minmove.apply(pose);

		TR << "CHI1\t" << pose.chi(1,2) << std::endl;
		TR << "CHI2\t" << pose.chi(2,2) << std::endl;
		TR << "CHI3\t" << pose.chi(3,2) << std::endl;
		core::Real const finalenergy((*sfxn)(pose));
		TR << "ENERGY\t" << finalenergy << std::endl;

		TS_ASSERT_DELTA( pose.chi(1,2), 29.90, 0.01 );
		TS_ASSERT_DELTA( pose.chi(2,2), 29.90, 0.01 );
		TS_ASSERT_DELTA( pose.chi(3,2), 29.90, 0.01 );
		TS_ASSERT_DELTA( finalenergy, 0.0, 0.01 );

		TR.flush();
	}

};



