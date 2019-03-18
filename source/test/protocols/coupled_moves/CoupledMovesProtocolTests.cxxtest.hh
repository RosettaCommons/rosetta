// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Amanda Loshbaugh (aloshbau@gmail.com)

#ifndef INCLUDED_protocols_coupled_moves_CoupledMovesProtocolTests_CXXTEST_HH
#define INCLUDED_protocols_coupled_moves_CoupledMovesProtocolTests_CXXTEST_HH

// Core headers
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelector.hh>
#include <core/pack/task/operation/ClashBasedRepackShell.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
//#include <core/conformation/Conformation.hh>
//#include <core/kinematics/FoldTree.hh>

// Protocol headers
#include <protocols/coupled_moves/CoupledMovesProtocol.hh>
#include <protocols/coupled_moves/CoupledMovesProtocolCreator.hh>

#include <basic/Tracer.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR("CoupledMovesProtocolTests_cxxtest");

using namespace std;
using core::import_pose::pose_from_file;
using core::import_pose::PDB_file;

class CoupledMovesProtocolTests : public CxxTest::TestSuite {

private:
	protocols::coupled_moves::CoupledMovesProtocolOP coupled_moves_protocol_;
	core::pack::task::TaskFactoryOP main_task_factory_;
	core::pose::PoseOP pose_ = utility::pointer::make_shared< core::pose::Pose >();

public:

	void setUp() {
		core_init_with_additional_options(
			"-extra_res_fa protocols/coupled_moves/input_files/NAD_from_1A80.params -resfile protocols/coupled_moves/input_files/1A80.resfile");
	}

	void test_apply() {
		coupled_moves_protocol_ = utility::pointer::make_shared< protocols::coupled_moves::CoupledMovesProtocol >();
		core::import_pose::pose_from_file( *pose_, "protocols/coupled_moves/input_files/1A80_with_NAD_truncated.pdb", core::import_pose::PDB_file );
		// Test that defaults are set up correctly
		// loop_size_
		core::Size loop_size = coupled_moves_protocol_->get_loop_size();
		TS_ASSERT( loop_size == 4 );
		// backbone_mover_
		std::string backbone_mover = coupled_moves_protocol_->get_backbone_mover();
		TS_ASSERT( backbone_mover == "backrub" );
		// repack_neighborhood_
		core::Real repack_neighborhood = coupled_moves_protocol_->get_repack_neighborhood();
		TS_ASSERT( repack_neighborhood == false );
		// ligand_mode_
		core::Real ligand_mode = coupled_moves_protocol_->get_ligand_mode();
		TS_ASSERT( ligand_mode == false );
		// number_ligands_
		core::Real number_ligands = coupled_moves_protocol_->get_number_ligands();
		TS_ASSERT( number_ligands == 1.0 );
		// ligand_weight_
		core::Real ligand_weight = coupled_moves_protocol_->get_ligand_weight();
		TS_ASSERT( ligand_weight == 1.0 );
		// ligand_prob_
		core::Real ligand_prob = coupled_moves_protocol_->get_ligand_prob();
		TS_ASSERT( ligand_prob == 0.1 );

		// Test that protocol produces both ligand and residue moves
		// and more residue than ligand moves
		coupled_moves_protocol_->set_ligand_mode( true );
		int ligand_moves = 0;
		int residue_moves = 0;
		core::Size iter = 50;
		for ( core::Size i = 1; i <= iter; ++i ) {
			core::Real const move_prob = numeric::random::uniform();
			coupled_moves_protocol_->setup_move_type( move_prob );
			std::string move_type = coupled_moves_protocol_->get_move_type();
			if ( move_type == "LIGAND" ) { ligand_moves = ( ligand_moves + 1 ); }
			else if ( move_type == "RESIDUE" ) { residue_moves = ( residue_moves + 1 ); }
		}
		TS_ASSERT( ligand_moves > 0 );
		TS_ASSERT( residue_moves > 0 );
		TS_ASSERT( residue_moves > ligand_moves );

		// run protocol
		coupled_moves_protocol_->set_ntrials( 10 );
		coupled_moves_protocol_->apply( *pose_ );

		// Test that protocol finds the correct move and design positions
		// based on the resfile. This also assumes that ClashBasedShellSelector is working
		utility::vector1<core::Size> move_positions = coupled_moves_protocol_->get_move_positions();
		utility::vector1<core::Size> design_positions = coupled_moves_protocol_->get_design_positions();
		TS_ASSERT_EQUALS(
			design_positions,
			utility::vector1<core::Size>({ 21, 41 })
		);
		TS_ASSERT_EQUALS(
			move_positions,
			utility::vector1<core::Size>({ 2, 14, 21, 22, 36, 41, 44, 46, 49, 64, 74, 96 })
		);

		// Test that ligand residue number has been correctly assigned
		utility::vector1<core::Size> ligand_resnums = coupled_moves_protocol_->get_ligand_resnums();
		TS_ASSERT_EQUALS(
			ligand_resnums,
			utility::vector1<core::Size>({ 96 })
		);

		// Test that protocol produces unique sequences
		std::map<std::string,core::Real> unique_sequences = coupled_moves_protocol_->get_unique_sequences();
		TR << unique_sequences << std::endl;
		core::Size n = unique_sequences.size();
		TS_ASSERT( n > 0 );
	} //test_apply

	///@brief check that improperly initialized CoupledMovesProtocol objects fail informatively
	void test_bad_init() {

		//test that ntrials is initialized
		protocols::coupled_moves::CoupledMovesProtocol cmp1;
		TS_ASSERT_EQUALS(cmp1.get_ntrials(), 1000);
		//we don't really care that this is 1000, mostly we care that it isn't an uninitialized variable, but 1000 is what it gets set to.

		return;
	}

}; //class CoupledMovesProtocolTests

#endif
