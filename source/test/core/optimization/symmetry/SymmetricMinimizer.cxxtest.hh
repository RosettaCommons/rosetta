// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/Minimizer.cxxtest.hh
/// @brief  test suite for Minimizer
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/symmetric_deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh> // RM: Needed for compilation on my platform.
#include <core/optimization/MinimizerOptions.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <core/pose/symmetry/util.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/Tracer.fwd.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.optimization.symmetry.SymmetricMinimizer.cxxtest");

using namespace core;

class SymmetricMinimizerTests : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	SymmetricMinimizerTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_c2_symm_score() {
		using namespace core;
		using namespace pose;
		using namespace conformation::symmetry;
		using namespace optimization;
		using namespace optimization::symmetry;

		core_init_with_additional_options( "-symmetry:symmetry_definition core/optimization/symmetry/c2.sdef -symmetry:initialize_rigid_body_dofs" );
		pose::Pose pose;
		///core::import_pose::pose_from_pdb( pose, "core/optimization/symmetry/c2_INPUT.pdb" );
		core::import_pose::pose_from_pdb( pose, "core/optimization/symmetry/c2_after_min.pdb" ); // the minimized structure gives better derivative agreement

		core::pose::symmetry::make_symmetric_pose( pose );

		scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();
		scoring::symmetry::SymmetricScoreFunctionOP sym_scorefxn = dynamic_cast< scoring::symmetry::SymmetricScoreFunction * > ( scorefxn() );

		//std::cout << "Symmetric score c2: " << (* scorefxn )( pose ) << std::endl;
		//std::cout << "Symmetric score before: " << (* scorefxn )( pose ) << std::endl;

		SymmetricConformation const & SymmConf ( dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
		assert( conformation::symmetry::is_symmetric( SymmConf ) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

		kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );

		//SymAtomTreeMinimizer minimizer;
		//MinimizerOptionsOP min_options = new MinimizerOptions( "dfpmin", 0.001, true, true, false );
		//minimizer.run( pose, *mm, *scorefxn, *min_options );

		SymmetricAtomDerivValidator sadv( pose, *sym_scorefxn, *mm );
		sadv.simple_deriv_check( true, 5e-4 );


		/*SymAtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options = new MinimizerOptions( "dfpmin", 0.001, true, true, false );
		minimizer.run( pose, *mm, *scorefxn, *min_options );

		pose.dump_pdb( "c2_after_min2.pdb" );*/

		//std::cout << "Symmetric score after: " << (* scorefxn )( pose ) << std::endl;
	}

	void test_c3_symm_score() {
		using namespace core;
		using namespace pose;
		using namespace scoring;
		using namespace conformation::symmetry;
		using namespace optimization;
		using namespace optimization::symmetry;

		core_init_with_additional_options( "-symmetry:symmetry_definition core/optimization/symmetry/c3.sdef -symmetry:initialize_rigid_body_dofs" );
		pose::Pose pose;
		///core::import_pose::pose_from_pdb( pose, "core/optimization/symmetry/c3_INPUT.pdb" );
		core::import_pose::pose_from_pdb( pose, "core/optimization/symmetry/c3_after_min.pdb" );
		core::pose::symmetry::make_symmetric_pose( pose );
		//pose.dump_pdb( "c3b_before_min.pdb" );

		scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();
		scoring::symmetry::SymmetricScoreFunctionOP sym_scorefxn = dynamic_cast< scoring::symmetry::SymmetricScoreFunction * > ( scorefxn() );

		//std::cout << "Symmetric score c3: " << (* scorefxn )( pose ) << std::endl;
		//std::cout << "Symmetric score before: " << (* scorefxn )( pose ) << std::endl;

		kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );

		SymmetricAtomDerivValidator sadv( pose, *sym_scorefxn, *mm );
		sadv.simple_deriv_check( true, 1e-2 );

		/// code below was used to take the c3_INPUT.pdb file and to minimize it into an interesting, collision free structure
		/*SymAtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options = new MinimizerOptions( "dfpmin", 0.00001, true, false, false );
		minimizer.run( pose, *mm, *scorefxn, *min_options );
		pose.dump_pdb( "c3_reminimized.pdb" );*/

		/*scorefxn->set_weight( fa_atr, 100 );
		scorefxn->set_weight( fa_rep, 0.02 );
		minimizer.run( pose, *mm, *scorefxn, *min_options );
		pose.dump_pdb( "c3b_mid_min1.pdb" );

		scorefxn->set_weight( fa_atr, 25 );
		scorefxn->set_weight( fa_rep, 0.05 );
		minimizer.run( pose, *mm, *scorefxn, *min_options );
		pose.dump_pdb( "c3b_mid_min2.pdb" );

		scorefxn->set_weight( fa_atr, 1 );
		scorefxn->set_weight( fa_rep, 0.1 );
		minimizer.run( pose, *mm, *scorefxn, *min_options );
		pose.dump_pdb( "c3b_mid_min3.pdb" );

		scorefxn->set_weight( fa_atr, 0.8 );
		scorefxn->set_weight( fa_rep, 0.2 );
		minimizer.run( pose, *mm, *scorefxn, *min_options );
		pose.dump_pdb( "c3b_mid_min4.pdb" );

		scorefxn->set_weight( fa_rep, 0.44 );
		minimizer.run( pose, *mm, *scorefxn, *min_options );

		pose.dump_pdb( "c3b_after_min.pdb" );*/

		//std::cout << "Symmetric score after: " << (* scorefxn )( pose ) << std::endl;
	}

	void test_fibril_symm_score() {
		using namespace core;
		using namespace pose;
		using namespace conformation::symmetry;
		using namespace optimization;
		using namespace optimization::symmetry;
		using namespace scoring;
		using namespace scoring::symmetry;

		core_init_with_additional_options( "-symmetry:symmetry_definition core/optimization/symmetry/fibril.symm -symmetry:initialize_rigid_body_dofs" );
		//core_init_with_additional_options( "-symmetry:symmetry_definition core/optimization/symmetry/fibril_min.symdef -symmetry:initialize_rigid_body_dofs" );

		//scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();
		//scoring::symmetry::SymmetricScoreFunctionOP sym_scorefxn = dynamic_cast< scoring::symmetry::SymmetricScoreFunction * > ( scorefxn() );

		scoring::EnergyMap weights;
		//1. etable
		weights[ fa_atr] = 0.8;
		weights[fa_rep ] = 0.44;
		weights[fa_sol ] = 0.6;
		verify_fibril_derivs_for_terms( weights, 3e-2 ); // lots of numerical noise from etable + etable derivs are wrong

		//2. hbonds
		weights.zero();
		weights[hbond_sr_bb ] = 1.17;
		weights[hbond_lr_bb ] = 1.17;
		weights[hbond_bb_sc ] = 1.10;
		weights[hbond_sc    ] = 1.10;
		verify_fibril_derivs_for_terms( weights, 1e-4 ); // Basically, the SER chi2 derivative problem means we have really bad derivatives.

		//3. pairE
		weights.zero();
		weights[fa_pair ] = 0.49;
		verify_fibril_derivs_for_terms( weights, 1e-5 );

		//4. dunE
		weights.zero();
		weights[fa_dun ] = 0.56;
		verify_fibril_derivs_for_terms( weights, 1e-5 );

		//5. proclose
		weights.zero();
		weights[pro_close ] = 1.0;
		verify_fibril_derivs_for_terms( weights, 1e-5 );

		//6. elec
		weights.zero();
		weights[ fa_elec ] = 0.6;
		verify_fibril_derivs_for_terms( weights, 1e-3 ); // lots of numerical noise from fa_elec
	}

	void verify_fibril_derivs_for_terms( scoring::EnergyMap const & in_weights, Real tolerance )
	{
		using namespace core;
		using namespace pose;
		using namespace conformation::symmetry;
		using namespace optimization;
		using namespace optimization::symmetry;
		using namespace scoring;
		using namespace scoring::symmetry;

		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "core/optimization/symmetry/3dg1_INPUT.pdb" );
		//core::import_pose::pose_from_pdb( pose, "core/optimization/symmetry/fibril_monomer.pdb" );
		core::pose::symmetry::make_symmetric_pose( pose );
		//pose.dump_pdb( "fibril_before_min.pdb" );

		SymmetricScoreFunctionOP sym_scorefxn = new SymmetricScoreFunction;
		scoring::ScoreFunctionOP scorefxn = sym_scorefxn;

		for ( Size ii = 1; ii <= (Size) end_of_score_type_enumeration; ++ii ) {
			if ( in_weights[ (ScoreType) ii ] != 0.0 ) {
				sym_scorefxn->set_weight( (ScoreType) ii, in_weights[ (ScoreType) ii ] );
			}
		}

		//std::cout << "Symmetric score fibril before: " << (* scorefxn )( pose ) << std::endl;
		//std::cout << "Symmetric score before: " << (* scorefxn )( pose ) << std::endl;

		SymmetricConformation const & SymmConf ( dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
		assert( conformation::symmetry::is_symmetric( SymmConf ) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

		kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );

		SymmetricAtomDerivValidator sadv( pose, *sym_scorefxn, *mm );
		//sadv.simple_deriv_check( true, 5e-2 );
		bool all_derivatives_correct = sadv.simple_deriv_check( true, tolerance );

		if ( ! all_derivatives_correct ) {
			//std::cout << "Derivative failure in the fibril symmetry for term(s): ";
			for ( Size ii = 1; ii <= (Size) end_of_score_type_enumeration; ++ii ) {
				if ( in_weights[ (ScoreType) ii ] != 0.0 ) {
					//std::cout  << " " << (ScoreType) ii;
				}
			}
			//std::cout << " given tolerance: " << tolerance << std::endl;
		}

		//SymAtomTreeMinimizer minimizer;
		//MinimizerOptionsOP min_options = new MinimizerOptions( "dfpmin", 0.00001, true, false, false );
		//minimizer.run( pose, *mm, *scorefxn, *min_options );
		//pose.dump_pdb( "fibril_minimized3.pdb" );

		//std::cout << "Symmetric score after: " << (* scorefxn )( pose ) << std::endl;

		/*kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
		SymAtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options = new MinimizerOptions( "dfpmin", 0.001, true, true, false );
		minimizer.run( pose, *mm, *scorefxn, *min_options );

		std::cout << "Symmetric score fibril after: " << (* scorefxn )( pose ) << std::endl;*/

	}

	void dont_test_helix_symm_score() {
		using namespace core;
		using namespace pose;
		using namespace conformation::symmetry;
		using namespace optimization;
		using namespace optimization::symmetry;

		core_init_with_additional_options( "-symmetry:symmetry_definition core/optimization/symmetry/helix.symm -symmetry:initialize_rigid_body_dofs" );
		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "core/optimization/symmetry/2bw9_INPUT.pdb" );
		core::pose::symmetry::make_symmetric_pose( pose );

		scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();

		//std::cout << "Symmetric score helix before: " << (* scorefxn )( pose ) << std::endl;

		//std::cout << "Symmetric score before: " << (* scorefxn )( pose ) << std::endl;

		SymmetricConformation const & SymmConf ( dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
		assert( conformation::symmetry::is_symmetric( SymmConf ) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

		kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );

		SymAtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options = new MinimizerOptions( "dfpmin", 0.001, true, true, false );
		minimizer.run( pose, *mm, *scorefxn, *min_options );

		//std::cout << "Symmetric score after: " << (* scorefxn )( pose ) << std::endl;

		/*kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
		SymAtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options = new MinimizerOptions( "dfpmin", 0.001, true, true, false );
		minimizer.run( pose, *mm, *scorefxn, *min_options );

		//std::cout << "Symmetric score helix after: " << (* scorefxn )( pose ) << std::endl;*/
	}



	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief simple test minimization
	void dont_test_simple_symmetric_min()
	{
		core_init_with_additional_options( "-symmetry:symmetry_definition core/scoring/symmetry/sym_def.dat" );

		using namespace optimization;
		using namespace optimization::symmetry;
		using pose::Pose;
		using id::AtomID;
		using id::DOF_ID;
		using id::PHI;
		using id::THETA;
		using id::D;

	//	pose::Pose start_pose(create_test_in_pdb_pose());
		pose::Pose start_pose;
		core::import_pose::pose_from_pdb( start_pose, "core/scoring/symmetry/test_in.pdb" );
		core::pose::symmetry::make_symmetric_pose( start_pose );

		kinematics::MoveMapOP mm = new kinematics::MoveMap;

		// setup moving dofs
		for ( int i=10; i<= 15; ++i ) {
			mm->set_bb ( i, true );
			mm->set_chi( i, true );
		}

		// setup the options
		scoring::ScoreFunctionOP scorefxn = new scoring::symmetry::SymmetricScoreFunction;

		SymAtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options = new MinimizerOptions( "linmin", 10.0, true, true, false );

	// non-core level code has no place in a "core" test!
	// furthermore, core does not (and should not) link the protocols lib!
	// move this test to protocols if you want to include/test code from protocols!!
	//	protocols::simple_moves::MinMover min_mover( mm, scorefxn, "linmin", 10.0, true
	//		/*use_nblist*/, true /*deriv_check*/, false /*no verbose-deriv-check, default*/ );

		{ // just fa_rama
			scorefxn->set_weight( scoring::rama, 1.0 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: rama" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // just fa_elec
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_elec, 0.5 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: fa_elec" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // just fa_dun
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_dun, 1.0 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: fa_dun" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // just fa_atr
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );

		}

		{ // just fa_atr, rep, sol
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr-rep-sol" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // fa_atr, rep, sol and fa_intra atr, rep, & sol
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );

			scorefxn->set_weight( scoring::fa_intra_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_intra_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_intra_sol, 0.65 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr-rep-sol and intra atr-rep-sol" << std::endl;
			TR << "start score: " << (*scorefxn)( pose ) << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
			pose.dump_pdb( "min_intrares.pdb" );
			TR << "end score: " << (*scorefxn)( pose ) << std::endl;
		}

		{  // p_aa_pp
			scorefxn->reset();
			scorefxn->set_weight( scoring::p_aa_pp, 0.29 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: p_aa_pp" << std::endl;
			TR << "start score: " << (*scorefxn)( pose ) << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
			pose.dump_pdb( "min_intrares.pdb" );
			TR << "end score: " << (*scorefxn)( pose ) << std::endl;
		}

		{ // fa_atr, rep, sol and p_aa_pp
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );
			scorefxn->set_weight( scoring::p_aa_pp, 0.29 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr-rep-sol and p_aa_pp" << std::endl;
			TR << "start score: " << (*scorefxn)( pose ) << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
			pose.dump_pdb( "min_intrares.pdb" );
			TR << "end score: " << (*scorefxn)( pose ) << std::endl;
		}


		{ // just fa_atr, rep, sol, rigid-body minimization
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );

			kinematics::MoveMapOP mm2 = new kinematics::MoveMap;
			mm2->set_jump( true );
			core::pose::symmetry::make_symmetric_movemap( start_pose, *mm2 );

			TR << "MINTEST: atr-rep-sol jumpmin" << std::endl;
			minimizer.run( start_pose, *mm2, *scorefxn, *min_options );

		}

		{ // just fa_pair
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_pair, 1.0 );

			TR << "MINTEST: fa_pair" << std::endl;
			minimizer.run( start_pose, *mm, *scorefxn, *min_options );
		}

		{ // just backbone hbonds
			scorefxn->reset();
			scorefxn->set_weight( scoring::hbond_lr_bb, 1.0 );
			scorefxn->set_weight( scoring::hbond_sr_bb, 1.0 );

			TR << "MINTEST: bb hbonds" << std::endl;
			minimizer.run( start_pose, *mm, *scorefxn, *min_options );
		}

		{ // all hbonds
			scorefxn->reset();
			scorefxn->set_weight( scoring::hbond_lr_bb, 1.0 );
			scorefxn->set_weight( scoring::hbond_sr_bb, 1.0 );
			scorefxn->set_weight( scoring::hbond_bb_sc, 1.0 );
			scorefxn->set_weight( scoring::hbond_sc, 1.0 );

			TR << "MINTEST: all hbonds" << std::endl;
			minimizer.run( start_pose, *mm, *scorefxn, *min_options );
		}

		return;

		/// EVERYTHING BELOW IS UNREACHABLE



		Pose pose;
		pose = start_pose;

		// set the moving dofs
		kinematics::MoveMapOP mm1 = new kinematics::MoveMap;
		kinematics::MoveMapOP mm2 = new kinematics::MoveMap;
		kinematics::MoveMapOP mm3 = new kinematics::MoveMap;
		kinematics::MoveMapOP mm4 = new kinematics::MoveMap;
		kinematics::MoveMapOP mm5 = new kinematics::MoveMap;
		// single backbone
		mm1->set_bb( 4, true );

		// all bb and chi
		mm2->set_bb( true );
		mm2->set_chi( true );

		// single dof
		mm3->set( DOF_ID( AtomID(1,4), PHI ), true );

		// everything!
		mm4->set( PHI, true );
		mm4->set( THETA, true );
		mm4->set( D, true );


		// everything! + neighborlist auto-update
		mm5->set( PHI, true );
		mm5->set( THETA, true );
		mm5->set( D, true );

		// setup scorefxn
		scorefxn->reset();
		scorefxn->set_weight( scoring::fa_atr, 0.80 );

		pose.dump_pdb( "before.pdb" );
		MinimizerOptionsOP min_options2 = new MinimizerOptions( "dfpmin", 0.001, true, true, false );
		minimizer.run( pose, *mm1, *scorefxn, *min_options2 );

		pose.dump_pdb( "after1.pdb" );

		exit(0);

		minimizer.run( pose, *mm2, *scorefxn, *min_options2 );
		pose.dump_pdb( "after2.pdb" );

		minimizer.run( pose, *mm3, *scorefxn, *min_options2 );
		pose.dump_pdb( "after3.pdb" );

		minimizer.run( pose, *mm4, *scorefxn, *min_options2 );
		pose.dump_pdb( "after4.pdb" );

		min_options2->nblist_auto_update( true );
		minimizer.run( pose, *mm5, *scorefxn, *min_options2 );
		pose.dump_pdb( "after5.pdb" );

	}


};
