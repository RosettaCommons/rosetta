// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/stepwise/modeler/rna//stepwise_modeler_rna_util_test.cxxtest.hh
/// @brief  test utils in stepwise/modeler/rna
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/modeler/rna/util.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("stepwise_modeler_rna_util_test");


class stepwise_modeler_rna_util_test : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-fasta protocols/stepwise/modeler/align/seq_with_AAmismatch.fasta -s protocols/stepwise/modeler/align/scaff_subset_stepwise_input_seq_1.pdb" );
	}

	void tearDown(){

	}



	void test_virtualize_free_rna_moieties(){

		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::stepwise::setup;
		using namespace protocols::stepwise::modeler::rna;
		using namespace utility;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		PoseOP pose_op = initialize_pose_and_other_poses_from_command_line( rsd_set );
		pose::Pose & pose = *pose_op;

		TS_ASSERT_EQUALS( pose.sequence(), "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTacugaggaucacccagu" );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "F[PHE:NtermProteinFull]ANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPI[ILE:CtermProteinFull]F[PHE:NtermProteinFull]ANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELT[THR:CtermProteinFull]a[RAD:LowerRNA:Virtual_Phosphate]cu[URA:Virtual_Phosphate]gaggaucacccag[RGU:Virtual_Phosphate]u[URA:UpperRNA]" );

		Size U_bulge = const_full_model_info( pose ).full_model_parameters()->conventional_to_full( 10, 'R' );
		Size A_bulge = const_full_model_info( pose ).full_model_parameters()->conventional_to_full(  6, 'R' );
		TS_ASSERT( !pose.residue_type( const_full_model_info( pose ).full_to_sub( U_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );
		TS_ASSERT( !pose.residue_type( const_full_model_info( pose ).full_to_sub( A_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );

		// check that virtualizing free RNA moieties bulges a bulged U.
		virtualize_free_rna_moieties( pose );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "F[PHE:NtermProteinFull]ANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPI[ILE:CtermProteinFull]F[PHE:NtermProteinFull]ANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELT[THR:CtermProteinFull]a[RAD:LowerRNA:Virtual_Phosphate]cu[URA:Virtual_Phosphate]gagga[RAD:3PrimePhos]u[URA:Virtual_RNA_Residue]c[RCY:Virtual_Phosphate]acccag[RGU:Virtual_Phosphate]u[URA:UpperRNA]" );
		TS_ASSERT(  pose.residue_type( const_full_model_info( pose ).full_to_sub( U_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );
		TS_ASSERT( !pose.residue_type( const_full_model_info( pose ).full_to_sub( A_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );


		// now just strip out RNA.
		utility::vector1< Size > rna_res;
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			if ( pose.residue_type( n ).is_RNA() ) rna_res.push_back( n );
		}
		pdbslice( pose, rna_res );

		TS_ASSERT_EQUALS( pose.sequence(), "acugaggaucacccagu" );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "a[RAD:LowerRNA:Virtual_Phosphate]cu[URA:Virtual_Phosphate]gagga[RAD:3PrimePhos]u[URA:Virtual_RNA_Residue]c[RCY:Virtual_Phosphate]acccag[RGU:Virtual_Phosphate]u[URA:UpperRNA]" );
		TS_ASSERT(  pose.residue_type( const_full_model_info( pose ).full_to_sub( U_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );
		TS_ASSERT( !pose.residue_type( const_full_model_info( pose ).full_to_sub( A_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );

		// check that virtualizing free RNA moieties virtualizes the bulged A.
		virtualize_free_rna_moieties( pose );

		// check that another A gets virtualized.
		// Actually a second one also gets virtualized but should not be. =( Let's not assert that it is bulged
		TS_ASSERT_EQUALS( pose.sequence(), "acugaggaucacccagu" );
		//  TS_ASSERT_EQUALS( pose.annotated_sequence(), "a[RAD:LowerRNA:Virtual_Phosphate]cu[URA:Virtual_Phosphate]ga[RAD:Virtual_RNA_Residue]g[RGU:Virtual_Phosphate]ga[RAD:3PrimePhos]u[URA:Virtual_RNA_Residue]c[RCY:Virtual_Phosphate]a[RAD:Virtual_RNA_Residue]c[RCY:Virtual_Phosphate]ccag[RGU:Virtual_Phosphate]u[URA:UpperRNA]" );
		TS_ASSERT(  pose.residue_type( const_full_model_info( pose ).full_to_sub( U_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );
		TS_ASSERT(  pose.residue_type( const_full_model_info( pose ).full_to_sub( A_bulge ) ).has_variant_type( VIRTUAL_RNA_RESIDUE ) );


	}



};



