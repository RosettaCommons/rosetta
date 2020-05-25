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
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.hh>
#include <core/import_pose/libraries/RNA_ChunkLibrary.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/pose/copydofs/CopyDofs.hh>


#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("RNA_DeNovoProtocolTest");

using namespace core;
using namespace protocols::rna::movers;
using namespace core::import_pose::options;
using namespace protocols::rna::denovo::movers;
using namespace core::import_pose::libraries;

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

};



