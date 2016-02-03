// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/denovo_design/BuildSheet.cxxtest.hh
/// @brief  test suite for devel::denovo_design::BuildSheet
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <devel/denovo_design/ParametricSheet.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>

/// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// C++ headers
#include <ObjexxFCL/FArray1D.hh>
#include <boost/lexical_cast.hpp>
#include <iostream>

static basic::Tracer TR("devel.denovo_design.BuildSheet.cxxtest");

// --------------- Test Class --------------- //
class BuildSheetTests : public CxxTest::TestSuite {
	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;

	// rotsets and ig to pass to packer
	core::pack::rotamer_set::RotamerSetsOP rotsets;
	core::pack::interaction_graph::AnnealableGraphBaseOP ig;

	// original pose read from disk
	core::pose::Pose input_pose;

public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();
		TR << "SETTING UP." << std::endl;
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if ( !residue_set.has_name("D2I") ) params_files.push_back("devel/denovo_design/D2I.params");
		residue_set.read_files_for_custom_residue_types(params_files);

		// initialize common filters/movers/scorefxns
		scorefxn = core::scoring::get_score_function( true );

		std::string const pdb_file( "devel/denovo_design/test_input.pdb" );
		core::import_pose::pose_from_file( input_pose, pdb_file , core::import_pose::PDB_file);

		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( input_pose ) );
		task->initialize_from_command_line();
		task->or_include_current(false);
		task->restrict_to_repacking();

		rotsets = core::pack::rotamer_set::RotamerSetsOP( new core::pack::rotamer_set::RotamerSets() );
		ig = core::pack::pack_rotamers_setup( input_pose, *scorefxn, task, rotsets );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// test overall acceptance criteria code
	void test_flat_sheet_grid() {
		// build a flat sheet -- result should be a planar, equally spaced grid
		typedef devel::denovo_design::ParametricSheet::StrandData StrandData;
		devel::denovo_design::ParametricSheet builder;
		builder.c_dist( 3.4 );
		builder.strand_dist( 4.5 );
		builder.c_coil( 0.0 );
		builder.h_coil( 0.0 );
		builder.twist( 0.0 );
		utility::vector1< StrandData > strand_data;
		for ( core::Size i=1; i<=4; ++i ) {
			strand_data.push_back( StrandData( "S" + boost::lexical_cast< std::string >( i ),
				6,
				0,
				"P" ) );
		}
		builder.strand_data( strand_data );
		// initialize points, ca_coords data
		builder.init_ca_coords();

		// there should be n+2 strands, each with length+4 entries
		// here all points should be -9999, -9999, -9999
		TS_ASSERT( builder.ca_coords_size() == 6 );
		for ( core::Size i=1; i<=builder.ca_coords_size(); ++i ) {
			TS_ASSERT( builder.ca_coords_size(i) == 10 );
			for ( core::Size j=1; j<=builder.ca_coords_size(i); ++j ) {
				TS_ASSERT_DELTA( builder.ca_coords(i,j).x(), -9999.9, 1e-6 );
				TS_ASSERT_DELTA( builder.ca_coords(i,j).y(), -9999.9, 1e-6 );
				TS_ASSERT_DELTA( builder.ca_coords(i,j).z(), -9999.9, 1e-6 );
			}
		}
		/*// strands should be build from the center outward
		for ( core::Size i=center_low; i>=1; --i ) {
		generate_strand( i, i+1, i+2 );
		}
		for ( core::Size i=center_high; i<=strands_to_build; ++i ) {
		generate_strand( i, i-1, i-2 );
		}*/

	}

};
