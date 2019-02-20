// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/PeptoidDesignTests.cxxtest.hh
/// @brief  Unit tests for designing with peptoids.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/pdb1rpb.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>

// Protocols Headers (for convenience)
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/aa_composition/AddCompositionConstraintMover.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("PeptoidDesignTests");


class PeptoidDesignTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}



	void test_design_001_601_602(){
		core::pose::Pose pose( pdb1rpb_pose() );
		std::stringstream selection_str;
		bool first (true);
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			if ( pose.residue_type(i).is_disulfide_bonded() ) continue;
			protocols::simple_moves::MutateResidue mutres( i, "GLY" );
			mutres.apply(pose); //Get rid of residue sidechains.
			if ( i == 1 ) continue;
			if ( !first ) {
				selection_str << ",";
			} else {
				first = false;
			}
			selection_str << i;
		}
		core::select::residue_selector::ResidueIndexSelectorOP index_selector( new core::select::residue_selector::ResidueIndexSelector );
		index_selector->set_index( selection_str.str() );
		core::select::residue_selector::NotResidueSelectorOP not_selector( new core::select::residue_selector::NotResidueSelector(index_selector) );

		using namespace core::pack::palette;
		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		core::pack::task::TaskFactoryOP task( new core::pack::task::TaskFactory );
		CustomBaseTypePackerPaletteOP pp = utility::pointer::make_shared< CustomBaseTypePackerPalette >();
		pp->add_base_residue_type( "001" );
		pp->add_base_residue_type( "601" );
		pp->add_base_residue_type( "602" );
		task->set_packer_palette( pp );
		core::pack::task::operation::ReadResfileOP resfile( new core::pack::task::operation::ReadResfile );
		resfile->filename( "core/pack/peptoid_design_1.resfile" );
		resfile->set_residue_selector( index_selector );
		task->push_back(resfile);

		core::pack::task::operation::RestrictToRepackingRLTOP restrict_to_repack( new core::pack::task::operation::RestrictToRepackingRLT );
		core::pack::task::operation::OperateOnResidueSubsetOP op_on_subset( new core::pack::task::operation::OperateOnResidueSubset( restrict_to_repack , not_selector ) );
		task->push_back(op_on_subset);

		protocols::minimization_packing::PackRotamersMover packer(sfxn, task->create_task_and_apply_taskoperations(pose));
		packer.apply(pose);

		TR << "AFTER PACKING:" << std::endl;
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			TR << "Position " << i << ": " << pose.residue_type(i).name() << std::endl;
			if ( i == 1 || pose.residue_type(i).is_disulfide_bonded() ) continue;
			TS_ASSERT( pose.residue_type(i).name3() == "001" || pose.residue_type(i).name3() == "601" || pose.residue_type(i).name3() == "602" );
		}
	}

	void test_design_hydrophobic_aa_001_601_602(){
		core::pose::Pose pose( pdb1rpb_pose() );
		std::stringstream selection_str;
		bool first (true);
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			if ( pose.residue_type(i).is_disulfide_bonded() ) continue;
			protocols::simple_moves::MutateResidue mutres( i, "GLY" );
			mutres.apply(pose); //Get rid of residue sidechains.
			if ( i == 1 ) continue;
			if ( !first ) {
				selection_str << ",";
			} else {
				first = false;
			}
			selection_str << i;
		}
		core::select::residue_selector::ResidueIndexSelectorOP index_selector( new core::select::residue_selector::ResidueIndexSelector );
		index_selector->set_index( selection_str.str() );
		core::select::residue_selector::NotResidueSelectorOP not_selector( new core::select::residue_selector::NotResidueSelector(index_selector) );

		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		sfxn->set_weight( core::scoring::aa_composition, 1.0 );
		protocols::aa_composition::AddCompositionConstraintMover add_comp;
		add_comp.create_constraint_from_file( "core/pack/peptoid_design_2.comp" );
		add_comp.apply(pose);

		using namespace core::pack::palette;
		core::pack::task::TaskFactoryOP task( new core::pack::task::TaskFactory );
		core::pack::task::operation::ReadResfileOP resfile( new core::pack::task::operation::ReadResfile );
		CustomBaseTypePackerPaletteOP pp = utility::pointer::make_shared< CustomBaseTypePackerPalette >();
		pp->add_base_residue_type( "001" );
		pp->add_base_residue_type( "601" );
		pp->add_base_residue_type( "602" );
		task->set_packer_palette( pp );
		resfile->filename( "core/pack/peptoid_design_2.resfile" );
		resfile->set_residue_selector( index_selector );
		task->push_back(resfile);

		core::pack::task::operation::RestrictToRepackingRLTOP restrict_to_repack( new core::pack::task::operation::RestrictToRepackingRLT );
		core::pack::task::operation::OperateOnResidueSubsetOP op_on_subset( new core::pack::task::operation::OperateOnResidueSubset( restrict_to_repack , not_selector ) );
		task->push_back(op_on_subset);

		protocols::minimization_packing::PackRotamersMover packer(sfxn, task->create_task_and_apply_taskoperations(pose));
		packer.apply(pose);

		TR << "AFTER PACKING:" << std::endl;
		//pose.dump_pdb( "QQQQQ.pdb" ); //DELETE ME
		bool at_least_one_aa(false), at_least_one_peptoid(false);
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			TR << "Position " << i << ": " << pose.residue_type(i).name() << std::endl;
			if ( i == 1 || pose.residue_type(i).is_disulfide_bonded() ) continue;
			if ( pose.residue_type(i).is_alpha_aa() ) {
				at_least_one_aa = true;
				core::chemical::AA const curaa( pose.residue_type(i).aa() );
				using namespace core::chemical;
				TS_ASSERT( curaa == aa_phe || curaa == aa_ala || curaa == aa_met || curaa == aa_ile || curaa == aa_tyr || curaa == aa_leu || curaa == aa_val || curaa == aa_trp );
			}
			if ( pose.residue_type(i).is_peptoid() ) {
				at_least_one_peptoid = true;
				TS_ASSERT( pose.residue_type(i).name3() == "001" || pose.residue_type(i).name3() == "601" || pose.residue_type(i).name3() == "602" );
			}
		}
		TS_ASSERT( at_least_one_aa );
		TS_ASSERT( at_least_one_peptoid );
	}

	/*void DISABLED_test_design_001_601_602_from_sarcosine(){
	core::pose::Pose pose( pdb1rpb_pose() );
	std::stringstream selection_str;
	bool first (true);
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
	if ( pose.residue_type(i).is_disulfide_bonded() ) continue;
	protocols::simple_moves::MutateResidue mutres( i, ( i == 1 ? "GLY" : "GLY:N_Methylation" ) );
	mutres.apply(pose); //Get rid of residue sidechains.
	if ( i == 1 ) continue;
	if ( !first ) {
	selection_str << ",";
	} else {
	first = false;
	}
	selection_str << i;
	}

	for ( core::Size i(2), imax(pose.total_residue()); i<=imax; ++i ) {
	if ( pose.residue_type(i).is_disulfide_bonded() ) continue;
	TS_ASSERT_EQUALS( pose.residue_type(i).name3(), "GLY" );
	if ( i != imax ) TS_ASSERT_EQUALS( pose.residue_type(i).name(), "GLY:N_Methylation" );
	TS_ASSERT_EQUALS( pose.residue_type(i).interchangeability_group(), "SAR" );
	}

	core::select::residue_selector::ResidueIndexSelectorOP index_selector( new core::select::residue_selector::ResidueIndexSelector );
	index_selector->set_index( selection_str.str() );
	core::select::residue_selector::NotResidueSelectorOP not_selector( new core::select::residue_selector::NotResidueSelector(index_selector) );

	core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
	core::pack::task::TaskFactoryOP task( new core::pack::task::TaskFactory );
	core::pack::task::operation::ReadResfileOP resfile( new core::pack::task::operation::ReadResfile );
	resfile->filename( "core/pack/peptoid_design_1.resfile" );
	resfile->set_residue_selector( index_selector );
	task->push_back(resfile);

	core::pack::task::operation::RestrictToRepackingRLTOP restrict_to_repack( new core::pack::task::operation::RestrictToRepackingRLT );
	core::pack::task::operation::OperateOnResidueSubsetOP op_on_subset( new core::pack::task::operation::OperateOnResidueSubset( restrict_to_repack , not_selector ) );
	task->push_back(op_on_subset);

	protocols::minimization_packing::PackRotamersMover packer(sfxn, task->create_task_and_apply_taskoperations(pose));
	packer.apply(pose);

	TR << "AFTER PACKING:" << std::endl;
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
	TR << "Position " << i << ": " << pose.residue_type(i).name() << std::endl;
	if ( i == 1 || pose.residue_type(i).is_disulfide_bonded() ) continue;
	TS_ASSERT( pose.residue_type(i).name3() == "001" || pose.residue_type(i).name3() == "601" || pose.residue_type(i).name3() == "602" );
	}
	}*/



};
