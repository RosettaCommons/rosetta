// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/rna/denovo//RNA_DeNovoProtocolTest.cxxtest.hh
/// @brief  RNA_DeNovoProtocol, wrapper for FARNA
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/sequence/util.hh>
#include <core/id/NamedAtomID.hh>

// Protocol Headers
#include <protocols/rna/denovo/RNA_DeNovoProtocol.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.hh>
#include <protocols/rna/movers/RNA_LoopCloser.hh>
#include <protocols/rna/denovo/options/RNA_DeNovoProtocolOptions.hh>
#include <protocols/rna/denovo/libraries/RNA_ChunkLibrary.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <core/pose/copydofs/CopyDofs.hh>


#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("RNA_DeNovoProtocolTest");

using namespace core;
using namespace protocols::rna::movers;
using namespace protocols::rna::denovo::options;
using namespace protocols::rna::denovo::movers;
using namespace protocols::rna::denovo::libraries;

class RNA_DeNovoProtocolTest : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-save_times -rna:denovo:cycles 1 -nstruct 1 -out:file:silent default.out -run:constant_seed" );
	}

	void tearDown(){

	}

	void test_save_times_option(){
		using namespace core::chemical;
		using namespace core::pose;
		using namespace protocols::rna::denovo;
		RNA_DeNovoProtocolOptionsOP options( new RNA_DeNovoProtocolOptions);
		options->initialize_from_command_line();
		RNA_DeNovoProtocol rna_de_novo_protocol( options );
		TS_ASSERT(  rna_de_novo_protocol.options()->save_times() );

		utility::file::file_delete( "default.out" );

		Pose pose;
		ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		make_pose_from_sequence( pose, "gggg", *rsd_set );
		TS_ASSERT_EQUALS( pose.size(), 4 );
		TS_ASSERT( pose.fold_tree().is_simple_tree() );

		rna_de_novo_protocol.apply( pose );
		TS_ASSERT_EQUALS( pose.size(), 4 );
		TS_ASSERT( pose.fold_tree().is_simple_tree() );
		TS_ASSERT( !rna_de_novo_protocol.rna_fragment_monte_carlo()->loop_modeling() );
		TS_ASSERT( hasPoseExtraScore( pose, "time" ) );
	}

	void test_sarcin_ricin_loop_fixed(){

		/////////////////////////////////////////////////////
		//
		//       _________________ Jump 1
		//      |                 |
		//   1  2 (3  4 5 x 6  7) 8  9  10
		//                               |
		//  19 18 (17 16x 15 14)  13 12 11
		//      |_________________| Jump 2
		//
		//  x mark possible cutpoints_closed
		//
		//  residues in parentheses will be built de novo
		//
		//  other residues come from starter PDB.
		//
		/////////////////////////////////////////////////////

		using namespace core;
		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::id;
		using namespace protocols::rna::denovo;
		using namespace protocols::toolbox;
		using namespace utility::tools;

		RNA_DeNovoProtocolOptionsOP options( new RNA_DeNovoProtocolOptions);
		options->initialize_from_command_line();
		// Sarcin-ricin loop modeling run
		options->set_rna_params_file("protocols/rna/farna_rebuild.params");
		options->set_chunk_pdb_files( make_vector1( "protocols/rna/srl_fixed_START1_1q9a_RNA.pdb") );
		options->set_input_res( utility::vector1< int >( std::get< 0 >( utility::get_resnum_and_chain( "1-2 8-13 18-19" ) ) ) );
		options->set_nstruct( 3 ); // by making > 1 models we actually test that atom_level_domain_map gets updated from model to model.
		RNA_DeNovoProtocol rna_de_novo_protocol( options );

		utility::file::file_delete( "default.out" );

		std::string const sequence = core::sequence::read_fasta_file_return_str( "protocols/rna/farna_rebuild.fasta" );
		Pose pose;
		ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		core::pose::make_pose_from_sequence( pose, sequence, *rsd_set );
		TS_ASSERT_EQUALS( pose.size(), 19 );
		TS_ASSERT_EQUALS( pose.sequence(), "ccuaguacgagaggaccgg" );
		TS_ASSERT( pose.fold_tree().is_simple_tree() );

		rna_de_novo_protocol.apply( pose );
		TS_ASSERT_EQUALS( pose.size(), 19 );
		TS_ASSERT( !pose.fold_tree().is_simple_tree() );
		TS_ASSERT_EQUALS( pose.fold_tree().num_jump(), 2 );
		TS_ASSERT( rna_de_novo_protocol.rna_fragment_monte_carlo()->loop_modeling() );

		AtomLevelDomainMap const & atom_level_domain_map = *(rna_de_novo_protocol.rna_fragment_monte_carlo()->atom_level_domain_map());
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 1 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 2 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 3 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 4 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 5 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 6 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 7 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 8 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 9 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",10 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",11 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",12 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",13 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",14 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",15 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",16 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",17 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",18 ), pose ), 1 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",19 ), pose ), 1 );

		atom_level_domain_map.show( TR );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " P  ", 1 ), pose ), core::pose::copydofs::FIXED_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " H5'", 1 ), pose ), core::pose::copydofs::FIXED_DOMAIN );

		bool found_cutpoint_closed_in_fixed_domain ( false );
		bool found_cutpoint_closed_in_moving_domain( false );
		utility::vector1< Size > moving_cutpoints_closed;
		for ( Size n = 1; n <= pose.size(); n++ ) {
			if ( pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) {
				TS_ASSERT( pose.fold_tree().is_cutpoint( n ) );
				TS_ASSERT( n < pose.size() );
				if ( atom_level_domain_map.get_domain( NamedAtomID( " C1'", n ), pose ) == 1 &&
						atom_level_domain_map.get_domain( NamedAtomID( " C1'", n ), pose ) == 1 ) {
					found_cutpoint_closed_in_fixed_domain = true;
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL1", n ), pose ), 1);
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL2", n ), pose ), 1);
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVU1", n+1 ), pose ), 1);
				} else {
					found_cutpoint_closed_in_moving_domain = true; // should be closing.
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL1", n ), pose ), 0);
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL2", n ), pose ), 0);
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVU1", n+1 ), pose ), 0);
					moving_cutpoints_closed.push_back( n );
				}
			}
		}
		TS_ASSERT( !found_cutpoint_closed_in_fixed_domain );
		TS_ASSERT( found_cutpoint_closed_in_moving_domain );
		TS_ASSERT_EQUALS( moving_cutpoints_closed.size(), 2 ); // one in each of the two loops.

		RNA_LoopCloser const & rna_loop_closer = *(rna_de_novo_protocol.rna_fragment_monte_carlo()->rna_loop_closer());
		utility::vector1< Size > loop_closer_cutpoints = rna_loop_closer.get_cutpoints_closed( pose );

		TS_ASSERT_EQUALS( loop_closer_cutpoints, moving_cutpoints_closed );
	}

	void test_sarcin_ricin_loop_free_with_base_pair_steps(){

		/////////////////////////////////////////////////////
		//
		//
		//   1 -2 x3 x4-5 x6 -7 x8-9  x 10
		//   |  |  |  |    |  |  |     /
		//  19 x18-17-16--15-14-13-12-11
		//
		//  x mark possible cutpoints_closed
		//
		/////////////////////////////////////////////////////

		using namespace core;
		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::id;
		using namespace protocols::rna::denovo;
		using namespace protocols::toolbox;
		using namespace utility::tools;

		RNA_DeNovoProtocolOptionsOP options( new RNA_DeNovoProtocolOptions);
		options->initialize_from_command_line();
		// Sarcin-ricin loop modeling run
		options->set_rna_params_file("protocols/rna/bps_test.params");
		options->set_nstruct( 2 ); // by making > 2 models we actually test that atom_level_domain_map gets updated from model to model.
		options->set_bps_moves( true );
		RNA_DeNovoProtocol rna_de_novo_protocol( options );

		utility::file::file_delete( "default.out" );

		std::string const sequence = core::sequence::read_fasta_file_return_str( "protocols/rna/farna_rebuild.fasta" );
		Pose pose;
		ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		core::pose::make_pose_from_sequence( pose, sequence, *rsd_set );
		TS_ASSERT_EQUALS( pose.size(), 19 );
		TS_ASSERT_EQUALS( pose.sequence(), "ccuaguacgagaggaccgg" );
		TS_ASSERT( pose.fold_tree().is_simple_tree() );

		rna_de_novo_protocol.apply( pose );
		TS_ASSERT_EQUALS( pose.size(), 19 );
		TS_ASSERT( !pose.fold_tree().is_simple_tree() );
		TS_ASSERT_EQUALS( pose.fold_tree().num_jump(), 7 );
		TS_ASSERT( !rna_de_novo_protocol.rna_fragment_monte_carlo()->loop_modeling() );

		AtomLevelDomainMap const & atom_level_domain_map = *(rna_de_novo_protocol.rna_fragment_monte_carlo()->atom_level_domain_map());
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 1 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 2 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 3 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 4 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 5 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 6 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 7 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 8 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 9 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",10 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",11 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",12 ), pose ), 0 );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",13 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",14 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",15 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",16 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",17 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",18 ), pose ), ROSETTA_LIBRARY_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",19 ), pose ), ROSETTA_LIBRARY_DOMAIN );

		atom_level_domain_map.show( TR );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " P  ", 1 ), pose ), core::pose::copydofs::FIXED_DOMAIN );
		TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " H5'", 1 ), pose ), core::pose::copydofs::FIXED_DOMAIN );

		bool found_cutpoint_closed_in_fixed_domain ( false );
		bool found_cutpoint_closed_in_moving_domain( false );
		utility::vector1< Size > moving_cutpoints_closed;
		for ( Size n = 1; n <= pose.size(); n++ ) {
			if ( pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) {
				TS_ASSERT( pose.fold_tree().is_cutpoint( n ) );
				TS_ASSERT( n < pose.size() );
				if ( atom_level_domain_map.get_domain( NamedAtomID( " C1'", n ), pose ) == 1 &&
						atom_level_domain_map.get_domain( NamedAtomID( " C1'", n ), pose ) == 1 ) {
					found_cutpoint_closed_in_fixed_domain = true;
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL1", n ), pose ), 1);
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL2", n ), pose ), 1);
					TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVU1", n+1 ), pose ), 1);
				} else {
					found_cutpoint_closed_in_moving_domain = true; // should be closing.
					// domain can also be ROSETTA_LIBRARY_DOMAIN
					// TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL1", n ), pose ), 0 );
					// TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVL2", n ), pose ), 0 );
					// TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( "OVU1", n+1 ), pose ), 0);
					moving_cutpoints_closed.push_back( n );
				}
			}
		}
		TS_ASSERT( !found_cutpoint_closed_in_fixed_domain );
		TS_ASSERT( found_cutpoint_closed_in_moving_domain );
		TS_ASSERT_EQUALS( moving_cutpoints_closed.size(), 7 );

		RNA_LoopCloser const & rna_loop_closer = *(rna_de_novo_protocol.rna_fragment_monte_carlo()->rna_loop_closer());
		utility::vector1< Size > loop_closer_cutpoints = rna_loop_closer.get_cutpoints_closed( pose );
		TS_ASSERT_EQUALS( loop_closer_cutpoints, moving_cutpoints_closed );
	}

};



