// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_pack_task_residue_selector_ClashBasedShellSelectorTests_CXXTEST_HH
#define INCLUDED_core_pack_task_residue_selector_ClashBasedShellSelectorTests_CXXTEST_HH

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>

// Core headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelector.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <set>

using namespace std;
using namespace core::pack::task::residue_selector;
using core::Size;
using core::Real;
using core::import_pose::pose_from_file;
using core::import_pose::PDB_file;
using core::pack::task::PackerTaskOP;
using core::pack::task::TaskFactory;
using core::pack::task::operation::OperateOnResidueSubset;
using core::pack::task::operation::OperateOnResidueSubsetOP;
using core::pack::task::operation::PreventRepackingRLT;
using core::pack::task::operation::ResLvlTaskOperationOP;
using core::pack::task::operation::RestrictToRepacking;
using core::pack::task::operation::RestrictToRepackingOP;
using core::pose::Pose;
using core::pose::make_pose_from_sequence;
using core::select::residue_selector::ResidueSelectorOP;
using core::select::residue_selector::ResidueIndexSelector;
using core::select::residue_selector::ResidueIndexSelectorOP;
using basic::datacache::DataMap;
using utility::vector1;


class ClashBasedShellSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options(
			"-extra_res_fa core/pack/task/residue_selector/NAD_from_1A80.params");
	}

	set<Size> make_shell(ResidueSelectorOP selector, Pose const & pose) {
		auto mask = selector->apply(pose);
		return resnums_from_bool_mask(mask);
	}

	set<Size> make_shell(set<Size> focus, Pose const & pose) {
		auto task = TaskFactory::create_packer_task( pose );
		for ( Size i=1; i <= pose.size(); i++ ) {
			if ( focus.count( i ) == 0 ) {
				task->nonconst_residue_task(i).restrict_to_repacking();
			}
		}
		auto selector = ClashBasedShellSelectorOP( new ClashBasedShellSelector( task ) );
		selector->set_include_focus(false);
		return make_shell(selector, pose);
	}

	set<Size> make_shell(string tag, DataMap & data, Pose const & pose) {
		auto selector = parse_selector_tag<ClashBasedShellSelector>(tag, data);
		return make_shell(selector, pose);
	}

	void print_clashes(
		utility::vector1< bool > clashing_positions,
		core::pose::Pose pose
	) {
		// print list of positions that clash for pymol
		std::cout << std::endl << "Clash-based repack shell:" << std::endl;
		core::Size num_positions=0;
		for ( core::Size i = 1; i <= clashing_positions.size(); ++i ) {
			if ( clashing_positions[i] ) {
				std::cout << pose.pdb_info()->number(i) << " "
					<< pose.pdb_info()->chain(i) << " NATAA" << std::endl;
				++num_positions;
			}
		}
		std::cout << std::endl << "For highlighting in pymol:" << std::endl;
		std::stringstream for_pymol;
		for ( core::Size i = 1; i <= clashing_positions.size(); ++i ) {
			if ( clashing_positions[i] ) {
				for_pymol << "chain " << pose.pdb_info()->chain(i) << " and resi " << pose.pdb_info()->number(i) << " or ";
			}
		}
		std::string for_pymol_string=for_pymol.str();
		std::cout << for_pymol_string.substr(0,for_pymol_string.size()-4) << std::endl;
		std::cout << std::endl << "Found " << num_positions << " positions that clash" << std::endl;

		// print list of positions that clash in format for unit test
		std::stringstream for_unit_test;
		for ( core::Size i = 1; i <= clashing_positions.size(); ++i ) {
			if ( clashing_positions[i] ) {
				for_unit_test << i << ",";
			}
		}
		for_unit_test << std::endl;
		std::string for_unit_test_string=for_unit_test.str();
		std::cout << for_unit_test_string;
	}

	void test_resnums_from_bool_mask() {
		// Empty mask:
		// This shouldn't cause any errors.
		TS_ASSERT_EQUALS(
			resnums_from_bool_mask( vector1<bool>({}) ),
			set<Size>({}) );

		// No true elements:
		TS_ASSERT_EQUALS(
			resnums_from_bool_mask({0}),
			set<Size>({}) );

		TS_ASSERT_EQUALS(
			resnums_from_bool_mask({0, 0}),
			set<Size>({}) );

		// One true element:
		// `resnums` should be 1-indexed.
		TS_ASSERT_EQUALS(
			resnums_from_bool_mask({1, 0}),
			set<Size>({1}) );

		TS_ASSERT_EQUALS(
			resnums_from_bool_mask({0, 1}),
			set<Size>({2}) );

		// All true elements:
		TS_ASSERT_EQUALS(
			resnums_from_bool_mask({1}),
			set<Size>({1}) );

		TS_ASSERT_EQUALS(
			resnums_from_bool_mask({1, 1}),
			set<Size>({1, 2}) );
	}

	void test_resnums_from_task() {
		Pose pose;
		PackerTaskOP task = TaskFactory::create_packer_task(pose);

		// Empty task
		// This shouldn't cause any errors.
		TS_ASSERT_EQUALS(
			resnums_from_task(task),
			set<Size>({}) );

		// Make a packer task for a pose with two residues.
		make_pose_from_sequence(pose, "AA", core::chemical::FA_STANDARD);
		task = TaskFactory::create_packer_task(pose);

		// Everything designable:
		// `resnums` should be 1-indexed.
		TS_ASSERT_EQUALS(
			resnums_from_task(task),
			set<Size>({1, 2}) );

		// Everything repackable, but not designable:
		task->restrict_to_repacking();
		TS_ASSERT_EQUALS(
			resnums_from_task(task),
			set<Size>({1, 2}) );

		// Nothing repackable:
		task->nonconst_residue_task(1).prevent_repacking();
		task->nonconst_residue_task(2).prevent_repacking();
		TS_ASSERT_EQUALS(
			resnums_from_task(task),
			set<Size>({}) );
	}

	void test_bool_mask_from_resnums() {
		Pose pose;

		// Empty pose:
		// This should cause any errors.
		TS_ASSERT_EQUALS(
			bool_mask_from_resnums(pose, {}),
			vector1<bool>({}) );

		// Make a pose with two residues.
		make_pose_from_sequence(pose, "AA", core::chemical::FA_STANDARD);

		// No residues:
		TS_ASSERT_EQUALS(
			bool_mask_from_resnums(pose, {}),
			vector1<bool>({0, 0}) );

		// One residue:
		TS_ASSERT_EQUALS(
			bool_mask_from_resnums(pose, {1}),
			vector1<bool>({1, 0}) );

		TS_ASSERT_EQUALS(
			bool_mask_from_resnums(pose, {2}),
			vector1<bool>({0, 1}) );

		// Two residues:
		TS_ASSERT_EQUALS(
			bool_mask_from_resnums(pose, {1, 2}),
			vector1<bool>({1, 1}) );
	}

	void test_bool_mask_from_task() {
		Pose pose;
		PackerTaskOP task = TaskFactory::create_packer_task(pose);

		// Empty task
		// This shouldn't cause any errors.
		TS_ASSERT_EQUALS(
			bool_mask_from_task(task),
			vector1<bool>({}) );

		// Make a packer task for a pose with two residues.
		make_pose_from_sequence(pose, "AA", core::chemical::FA_STANDARD);
		task = TaskFactory::create_packer_task(pose);

		// Everything designable:
		TS_ASSERT_EQUALS(
			bool_mask_from_task(task),
			vector1<bool>({1, 1}) );

		// Everything repackable, but not designable:
		task->restrict_to_repacking();
		TS_ASSERT_EQUALS(
			bool_mask_from_task(task),
			vector1<bool>({1, 1}) );

		// Nothing repackable:
		task->nonconst_residue_task(1).prevent_repacking();
		task->nonconst_residue_task(2).prevent_repacking();
		TS_ASSERT_EQUALS(
			bool_mask_from_task(task),
			vector1<bool>({0, 0}) );
	}

	void test_task_from_resnums() {
		Pose pose;
		PackerTaskOP task;

		// Empty pose:
		// This shouldn't cause any errors.
		task = task_from_resnums(pose, {});
		TS_ASSERT_EQUALS( task->total_residue(), 0 );

		// Make a pose with two residues.
		make_pose_from_sequence(pose, "AA", core::chemical::FA_STANDARD);

		// No residues:
		task = task_from_resnums(pose, {});
		TS_ASSERT_EQUALS( task->being_packed(1), false );
		TS_ASSERT_EQUALS( task->being_packed(2), false );

		// One residue:
		task = task_from_resnums(pose, {1});
		TS_ASSERT_EQUALS( task->being_packed(1), true );
		TS_ASSERT_EQUALS( task->being_packed(2), false );

		task = task_from_resnums(pose, {2});
		TS_ASSERT_EQUALS( task->being_packed(1), false );
		TS_ASSERT_EQUALS( task->being_packed(2), true );

		// Two residues:
		task = task_from_resnums(pose, {1, 2});
		TS_ASSERT_EQUALS( task->being_packed(1), true );
		TS_ASSERT_EQUALS( task->being_packed(2), true );
	}

	void test_task_from_bool_mask() {
		Pose pose;
		PackerTaskOP task;

		// Empty pose:
		// This shouldn't cause any errors.
		task = task_from_bool_mask(pose, vector1<bool>({}) );
		TS_ASSERT_EQUALS( task->total_residue(), 0 );

		// It's an error for the pose and the mask to have different lengths, but we
		// can't test for it because runtime_assert calls std::exit().

		// Make a pose with two residues.
		make_pose_from_sequence(pose, "AA", core::chemical::FA_STANDARD);

		// No true elements:
		task = task_from_bool_mask(pose, {0, 0});
		TS_ASSERT_EQUALS( task->being_packed(1), false );
		TS_ASSERT_EQUALS( task->being_packed(2), false );

		// One true element:
		task = task_from_bool_mask(pose, {1, 0});
		TS_ASSERT_EQUALS( task->being_packed(1), true );
		TS_ASSERT_EQUALS( task->being_packed(2), false );

		task = task_from_bool_mask(pose, {0, 1});
		TS_ASSERT_EQUALS( task->being_packed(1), false );
		TS_ASSERT_EQUALS( task->being_packed(2), true );

		// Two residues:
		task = task_from_bool_mask(pose, {1, 1});
		TS_ASSERT_EQUALS( task->being_packed(1), true );
		TS_ASSERT_EQUALS( task->being_packed(2), true );
	}

	void test_is_sc_xx_clash() {
		Pose no_clash, sc_sc_clash, sc_bb_clash;
		pose_from_file(no_clash, "core/pack/task/residue_selector/no_clash.pdb", PDB_file);
		pose_from_file(sc_sc_clash, "core/pack/task/residue_selector/sc_sc_clash.pdb", PDB_file);
		pose_from_file(sc_bb_clash, "core/pack/task/residue_selector/sc_bb_clash.pdb", PDB_file);

		// No clashes:
		TS_ASSERT_EQUALS(
			is_sc_sc_clash(no_clash.residue(1), no_clash.residue(2), 0.5),
			false);
		TS_ASSERT_EQUALS(
			is_sc_sc_clash(no_clash.residue(2), no_clash.residue(1), 0.5),
			false);
		TS_ASSERT_EQUALS(
			is_sc_bb_clash(no_clash.residue(1), no_clash.residue(2), 0.5),
			false);
		TS_ASSERT_EQUALS(
			is_sc_bb_clash(no_clash.residue(2), no_clash.residue(1), 0.5),
			false);

		// Just sidechain clashes:
		TS_ASSERT_EQUALS(
			is_sc_sc_clash(sc_sc_clash.residue(1), sc_sc_clash.residue(2), 0.5),
			true);
		TS_ASSERT_EQUALS(
			is_sc_sc_clash(sc_sc_clash.residue(2), sc_sc_clash.residue(1), 0.5),
			true);
		TS_ASSERT_EQUALS(
			is_sc_bb_clash(sc_sc_clash.residue(1), sc_sc_clash.residue(2), 0.5),
			false);
		TS_ASSERT_EQUALS(
			is_sc_bb_clash(sc_sc_clash.residue(2), sc_sc_clash.residue(1), 0.5),
			false);

		// Sidechain and backbone clashes:
		TS_ASSERT_EQUALS(
			is_sc_sc_clash(sc_bb_clash.residue(1), sc_bb_clash.residue(2), 0.5),
			true);
		TS_ASSERT_EQUALS(
			is_sc_sc_clash(sc_bb_clash.residue(2), sc_bb_clash.residue(1), 0.5),
			true);
		TS_ASSERT_EQUALS(
			is_sc_bb_clash(sc_bb_clash.residue(1), sc_bb_clash.residue(2), 0.5),
			true);
		TS_ASSERT_EQUALS(
			is_sc_bb_clash(sc_bb_clash.residue(2), sc_bb_clash.residue(1), 0.5),
			false);
	}

	void test_add_clashes_to_shell() {
		set<Size> focus, shell;
		Pose no_clash, sc_sc_clash, sc_bb_clash;
		pose_from_file(no_clash, "core/pack/task/residue_selector/no_clash.pdb", PDB_file);
		pose_from_file(sc_sc_clash, "core/pack/task/residue_selector/sc_sc_clash.pdb", PDB_file);
		pose_from_file(sc_bb_clash, "core/pack/task/residue_selector/sc_bb_clash.pdb", PDB_file);

		// No clashes.
		shell.clear();
		add_clashes_to_shell(no_clash, no_clash.residue(1), 0.5, focus, shell);
		TS_ASSERT_EQUALS( shell, set<Size>({}) );

		// No clashes.
		shell.clear();
		add_clashes_to_shell(no_clash, no_clash.residue(1), 0.5, focus, shell);
		TS_ASSERT_EQUALS( shell, set<Size>({}) );

		// Just sidechain clashes:
		shell.clear();
		add_clashes_to_shell(sc_sc_clash, sc_sc_clash.residue(1), 0.5, focus, shell);
		TS_ASSERT_EQUALS( shell, set<Size>({2}) );

		shell.clear();
		add_clashes_to_shell(sc_sc_clash, sc_sc_clash.residue(2), 0.5, focus, shell);
		TS_ASSERT_EQUALS( shell, set<Size>({1}) );

		// Sidechain and backbone clashes:
		shell.clear();
		add_clashes_to_shell(sc_bb_clash, sc_bb_clash.residue(1), 0.5, focus, shell);
		TS_ASSERT_EQUALS( shell, set<Size>({}) );

		// Don't count residues that are already part of the focus.
		shell.clear(); focus = {2};
		add_clashes_to_shell(sc_sc_clash, sc_sc_clash.residue(1), 0.5, focus, shell);
		TS_ASSERT_EQUALS( shell, set<Size>({}) );
	}

	void test_selections() {
		Pose pose1, pose2;
		pose_from_file(pose1, "core/pack/task/residue_selector/1A80_with_NAD.pdb", PDB_file);
		pose_from_file(pose2, "core/pack/task/residue_selector/2NXX.pdb", PDB_file);

		TS_ASSERT_EQUALS(
			make_shell({5}, pose1),
			set<Size>({ 3,7,15,182,227 }) );
		TS_ASSERT_EQUALS(
			make_shell({104}, pose1),
			set<Size>({ 73,75,89,93,101,123,126,132,135 }) );
		TS_ASSERT_EQUALS(
			make_shell({236}, pose1),
			set<Size>({ 240 }) );

		// Test if CBRSS can find
		// - clashing positions in multiple chains
		// - clashes with residue/ligands defined by -extra_res_fa ligand.params
		TS_ASSERT_EQUALS(
			make_shell({108}, pose1),
			set<Size>({ 76,77,107,109,110,139,273,276,278 }) );

		// Test if CBRSS can find clashes when
		// - the input residue is in chain A (the 1st chain) and should interact with residues in both chain A and chain E (the 2nd chain)
		TS_ASSERT_EQUALS(
			make_shell({193}, pose2),
			set<Size>({ 74,197,423,424,426,427 }) );

		// Test if CBRSS can find clashes when
		// - the input residue is in chain E (the 2nd chain) and should interact with residues in both chain A and chain E
		TS_ASSERT_EQUALS(
			make_shell({424}, pose2),
			set<Size>({ 144,146,160,193,352,353,356,357,421,427,428 }) );

		// Test if CBRSS returns no clashes when you give it no input
		TS_ASSERT_EQUALS(
			make_shell(set<Size>({}), pose2),
			set<Size>({}) );

		// KBK: I removed tests involving residues that are not in the pose, because
		// this is now an error.

		// Test if CBRSS works correctly if you give it multiple positions
		TS_ASSERT_EQUALS(
			make_shell({108,104}, pose1),
			set<Size>({ 73,75,76,77,89,93,101,107,109,110,123,126,132,135,139,273,276,278 }) );

		// Test if CBRSS works correctly if you give it the entire pose
		set<Size> whole_pose;
		for ( Size i = 1; i <= pose2.size(); ++i ) whole_pose.insert( i );
		TS_ASSERT_EQUALS(
			make_shell(whole_pose, pose2),
			set<Size>({}) );
	}

	void test_parse_my_tag() {
		Pose pose;
		pose_from_file(pose, "core/pack/task/residue_selector/1A80_with_NAD.pdb", PDB_file);

		// I had to take a closer look at the shell for residue 104 because it
		// excludes some positions I initially expected it to include, based on
		// looking at rotamers generated by pymol.  It turns out that the shell is
		// correct.  Here are the positions I initially thought would be included and
		// the reasons why they're not:
		//
		// 101: The rotamers that clash with this position reach all the way through
		// the sidechain and clash with the backbone.
		//
		// 89: The rotamers that clash with this position also clash with the
		// backbone of residue 101.
		//
		// 132+135: Pymol has a rotamer pointing towards these two residues, but
		// rosetta doesn't.

		auto sele = ResidueIndexSelectorOP( new ResidueIndexSelector );
		sele->set_index("104");

		auto taskop = OperateOnResidueSubsetOP( new OperateOnResidueSubset );
		taskop->op( ResLvlTaskOperationOP( new PreventRepackingRLT ) );
		taskop->selector( sele );
		taskop->flip_subset(true);

		basic::datacache::DataMap data;
		data.add("ResidueSelector", "sel_104", sele);
		data.add("task_operations", "op_104", taskop);
		data.add("task_operations", "no_design", RestrictToRepackingOP( new RestrictToRepacking ) );

		// Empty selector; no errors.
		TS_ASSERT_EQUALS(
			make_shell("<ClashBasedShell/>", data, pose),
			set<Size>({}) );

		// String residue list.
		TS_ASSERT_EQUALS(
			make_shell(R"(<ClashBasedShell resnums="104"/>)", data, pose),
			set<Size>({ 75, 104, 123 }) );

		// Residue selectors
		TS_ASSERT_EQUALS(
			make_shell(R"(<ClashBasedShell residue_selector="sel_104"/>)", data, pose),
			set<Size>({ 75, 104, 123 }) );

		TS_ASSERT_EQUALS(
			make_shell(
			R"(<ClashBasedShell>)"
			R"(  <Index resnums="104"/>)"
			R"(</ClashBasedShell>)", data, pose),
			set<Size>({ 75, 104, 123 }) );

		// Task operations
		// More positions selected because the task operation allows design.
		TS_ASSERT_EQUALS(
			make_shell(R"(<ClashBasedShell task_operations="op_104"/>)", data, pose),
			set<Size>({ 73, 75, 89, 93, 101, 104, 123, 126, 132, 135 }) );  // More positions because design is allowed.

		// Without design, and with focusing on residues that can only be repacked,
		// we get the same shell as before.
		TS_ASSERT_EQUALS(
			make_shell(R"(<ClashBasedShell task_operations="op_104,no_design" focus_on_designable="no"/>)", data, pose),
			set<Size>({ 75, 104, 123 }) );

		// Exclude residue 104 itself.
		TS_ASSERT_EQUALS(
			make_shell(R"(<ClashBasedShell resnums="104" include_focus="no"/>)", data, pose),
			set<Size>({ 75, 123 }) );

		// Two shells
		TS_ASSERT_EQUALS(
			make_shell(R"(<ClashBasedShell resnums="104" num_shells="2"/>)", data, pose),
			set<Size>({ 75, 104, 106, 120, 123, 135 }) );
	}


};

#endif
