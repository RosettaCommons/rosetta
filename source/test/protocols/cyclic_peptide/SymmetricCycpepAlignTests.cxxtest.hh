// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/SymmetricCycpepAlignTests.cxxtest.hh
/// @brief  Unit tests for the SymmetricCycpepAlign mover, which aligns a quasi-symmetric
/// peptide's symmetry axis with the z-axis and centers it on the origin.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/cyclic_peptide/SymmetricCycpepAlign.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

static THREAD_LOCAL basic::Tracer TR("SymmetricCycpepAlignTests");


class SymmetricCycpepAlignTests : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init_with_additional_options("-write_all_connect_info -extra_res_fa crosslinker/1.3.5_trisbromomethylbenzene.params sidechain_conjugation/CYX.params");

		core::pose::PoseOP pose_asymm( core::import_pose::pose_from_file( "protocols/cyclic_peptide/asymm.pdb", false, core::import_pose::PDB_file) );
		protocols::cyclic_peptide::DeclareBond bond0;
		bond0.set( 8, "C", 1, "N", false );
		bond0.apply( *pose_asymm );
		pose_asymm_ = pose_asymm;

		core::pose::PoseOP pose_c4m( core::import_pose::pose_from_file( "protocols/cyclic_peptide/c4m_symm.pdb", false, core::import_pose::PDB_file) );
		protocols::cyclic_peptide::DeclareBond bond1;
		bond1.set( 20, "C", 1, "N", false );
		bond1.apply( *pose_c4m );
		pose_c4m_ = pose_c4m;

		core::pose::PoseOP pose_c3_tbmb( core::import_pose::pose_from_file( "protocols/cyclic_peptide/c3_symm_tbmb.pdb", false, core::import_pose::PDB_file) );
		protocols::simple_moves::MutateResidue mut1(4, "CYX"), mut2(10, "CYX"), mut3(16, "CYX");
		mut1.apply( *pose_c3_tbmb );
		mut2.apply( *pose_c3_tbmb );
		mut3.apply( *pose_c3_tbmb );
		protocols::cyclic_peptide::DeclareBond bond2, bond3, bond4, bond5;
		bond2.set( 18, "C", 1, "N", false );
		bond3.set( 10, "SG", 19, "CM1", false );
		bond4.set( 16, "SG", 19, "CM2", false );
		bond5.set( 4, "SG", 19, "CM3", false );
		bond2.apply( *pose_c3_tbmb );
		bond3.apply( *pose_c3_tbmb );
		bond4.apply( *pose_c3_tbmb );
		bond5.apply( *pose_c3_tbmb );
		pose_c3_tbmb_ = pose_c3_tbmb;

	}

	void tearDown(){

	}


	/// @brief Tests whether we can auto-detect c4/m quasi-symmetry in a quasi-symmetric c4/m peptide pose.
	///
	void test_detect_c4m_symmetry() {
		TR << "Running SymmetricCycpepAlignTests::test_detect_c4m_symmetry() to check whether we can auto-detect c4/m symmetry in a quasi-symmetric peptide." << std::endl;
		core::pose::PoseOP pose( pose_c4m_->clone() );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover;
		mover.set_auto_detect_symmetry(true);
		mover.set_angle_threshold(15);

		mover.apply( *pose );

		TS_ASSERT_EQUALS( mover.last_symmetry_repeats(), 4);
		TS_ASSERT_EQUALS( mover.last_symmetry_mirror(), true);
	}

	/// @brief Tests whether we can center a c4/m-symmetric peptide on the origin.
	///
	void test_align_c4m_to_origin(){
		TR << "Running SymmetricCycpepAlignTests::test_align_c4m_to_origin() to check whether a c4/m quasi-symmetric peptide is properly aligned to the origin by the mover." << std::endl;
		core::pose::PoseOP pose( pose_c4m_->clone() );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover;
		mover.set_auto_detect_symmetry(true);
		mover.set_angle_threshold(15);

		mover.apply( *pose );
		numeric::xyzVector <core::Real> avgpos;
		avgpos.x(0); avgpos.y(0); avgpos.z(0);
		core::Size counter(0);
		for ( core::Size ir=1, irmax=pose->total_residue(); ir<=irmax; ++ir ) {
			if ( !pose->residue_type(ir).is_protein() ) continue;
			for ( core::Size ia=1; ia<5; ++ia ) {
				avgpos += pose->residue(ir).xyz(ia);
				++counter;
			}
		}

		avgpos /= static_cast< core::Real >( counter );

		TS_ASSERT_DELTA( avgpos.x(), 0, 1e-6 );
		TS_ASSERT_DELTA( avgpos.y(), 0, 1e-6 );
		TS_ASSERT_DELTA( avgpos.z(), 0, 1e-6 );

	}

	void test_align_c4m_to_zaxis(){
		TR << "Running SymmetricCycpepAlignTests::test_align_c4m_to_origin() to check whether a c4/m quasi-symmetric peptide is properly aligned to the z-axis by the mover." << std::endl;
		core::pose::PoseOP pose( pose_c4m_->clone() );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover;
		mover.set_auto_detect_symmetry(true);
		mover.set_angle_threshold(15);

		mover.apply( *pose );

		numeric::xyzVector< core::Real > v1( pose->residue(1).xyz(2) - pose->residue(11).xyz(2) );
		numeric::xyzVector< core::Real > v2( pose->residue(4).xyz(2) - pose->residue(14).xyz(2) );
		numeric::xyzVector< core::Real > crosspdt( v1.cross(v2).normalize() );
		numeric::xyzVector< core::Real > zaxis( 0, 0, 1 );

		core::Real const dotpdt( zaxis.dot( crosspdt ) );
		TR << "Dot product of peptide axis with z-axis is " << dotpdt << ".  (Expect 1.0.)" << std::endl;

		TS_ASSERT_DELTA( dotpdt, 1.0, 1e-4 );

	}

	void test_trim_c4m_to_repeat_two(){
		TR << "Running SymmetricCycpepAlignTests::test_trim_c4m_to_repeat_one() to check whether a c4/m quasi-symmetric peptide is properly trimmed after alignment." << std::endl;
		core::pose::PoseOP pose( pose_c4m_->clone() ), pose2( pose_c4m_->clone() );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover;
		mover.set_auto_detect_symmetry(true);
		mover.set_angle_threshold(15);
		mover.apply( *pose );

		mover.set_trim_info( true, 2 );
		mover.apply( *pose2 );

		TS_ASSERT_EQUALS( pose2->total_residue(), 5);

		for ( core::Size ir=1; ir<=5; ++ir ) {
			TS_ASSERT_DELTA( pose2->residue(ir).xyz("CA"), pose->residue(ir+5).xyz("CA"), 1e-6 );
		}
	}

	void test_detect_asymm(){
		TR << "Running SymmetricCycpepAlignTests::test_detect_asymm() to confirm that the mover can detect asymmetry or the wrong symmetry." << std::endl;
		core::pose::PoseOP pose( pose_asymm_->clone() ), pose2( pose_c3_tbmb_->clone() );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover;
		mover.set_auto_detect_symmetry(true);
		mover.set_angle_threshold(35);
		mover.apply( *pose );
		TS_ASSERT( mover.get_last_move_status() == protocols::moves::FAIL_DO_NOT_RETRY );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover2;
		mover2.set_symmetry(6, true);
		mover2.set_angle_threshold(35);
		mover2.apply( *pose2 );
		TS_ASSERT( mover2.get_last_move_status() == protocols::moves::FAIL_DO_NOT_RETRY );

	}

	void test_detect_c3_symmetry_tbmb(){
		TR << "Running SymmetricCycpepAlignTests::test_detect_c4m_symmetry() to check whether we can auto-detect c3 symmetry in a quasi-symmetric peptide with TBMB." << std::endl;
		core::pose::PoseOP pose( pose_c3_tbmb_->clone() );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover;
		mover.set_auto_detect_symmetry(true);
		mover.set_angle_threshold(35);

		mover.apply( *pose );

		TS_ASSERT_EQUALS( mover.last_symmetry_repeats(), 3);
		TS_ASSERT_EQUALS( mover.last_symmetry_mirror(), false);
	}

	void test_trim_c3_tbmb_to_repeat_two(){
		TR << "Running SymmetricCycpepAlignTests::test_trim_c4m_to_repeat_two() to check whether a c3 quasi-symmetric peptide with TBMB is properly trimmed after alignment." << std::endl;
		core::pose::PoseOP pose( pose_c3_tbmb_->clone() ), pose2( pose_c3_tbmb_->clone() );

		protocols::cyclic_peptide::SymmetricCycpepAlign mover;
		mover.set_auto_detect_symmetry(true);
		mover.set_angle_threshold(35);
		mover.apply( *pose );

		mover.set_trim_info( true, 2 );
		mover.apply( *pose2 );

		TS_ASSERT_EQUALS( pose2->total_residue(), 6);

		for ( core::Size ir=1; ir<=6; ++ir ) {
			TS_ASSERT_DELTA( pose2->residue(ir).xyz("CA"), pose->residue(ir+6).xyz("CA"), 1e-6 );
		}
	}

private:

	/// @brief A test pose, asymmetric.
	/// @details Clone this to manipulate.
	core::pose::PoseCOP pose_asymm_;

	/// @brief A test pose, c4/m symmetric.
	/// @details Clone this to manipulate.
	core::pose::PoseCOP pose_c4m_;

	/// @brief A test pose, c3 symmetric with tbmb.
	/// @details Clone this to manipulate.
	core::pose::PoseCOP pose_c3_tbmb_;


};



