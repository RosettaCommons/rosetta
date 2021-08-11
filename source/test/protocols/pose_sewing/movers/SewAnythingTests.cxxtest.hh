// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pose_sewing/movers/SewAnythingTests.cxxtest.hh
/// @brief  tests the full SewAnything pipeline to avoid needing an external segment file
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>

// Project Headers
#include <protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.hh>
#include <protocols/pose_sewing/movers/PickRandomSegmentMover.hh>
#include <protocols/pose_sewing/movers/AddFlankingVirtualResiduesMover.hh>
#include <protocols/pose_sewing/movers/SewAnythingAddMover.hh>
#include <core/select/residue_selector/SecondaryStructureSelector.hh>
#include <core/select/residue_selector/VirtualResidueSelector.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("SewAnythingTests");


class SewAnythingTests : public CxxTest::TestSuite {
	utility::vector1< protocols::filters::FilterOP > filters_;
	utility::vector1< core::simple_metrics::SimpleMetricOP > metrics_;
	//Define Variables

public:

	void setUp() {
		core_init();
		core_init_with_additional_options( "-corrections::beta_nov16 true" );
	}
	void test_empty(){
		protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVectorOP seg_vectorOP(new protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVector);
		seg_vectorOP->clear_segment_envelopes();
		seg_vectorOP->add_segment_envelopes_from_string("H 4 100 L 1 100 H 4 100");
		core::pose::Pose test_3_bundle = create_3_bundle_pose();
		std::string segment_file_path = "test_segments";
		utility::file::file_delete( "test_segments/segment.txt" );
		if ( ! utility::file::is_directory(segment_file_path) ) {
			utility::file::create_directory_recursive(segment_file_path);
		}

		seg_vectorOP->simultaneously_populate_and_write_segment_file_pdb_using_stored_motifs ( test_3_bundle, segment_file_path,0,1,filters_,metrics_,false);

		protocols::pose_sewing::movers::PickRandomSegmentMoverOP pick_mover = protocols::pose_sewing::movers::PickRandomSegmentMoverOP(new protocols::pose_sewing::movers::PickRandomSegmentMover());
		utility::vector1<std::string> paths;
		paths.push_back(segment_file_path);
		pick_mover->set_segment_file_paths(paths);
		core::Size orig_size = test_3_bundle.size();
		pick_mover->apply(test_3_bundle);
		TS_ASSERT(test_3_bundle.size() < orig_size);
		orig_size = test_3_bundle.size();

		core::select::residue_selector::SecondaryStructureSelectorOP in_selector( new core::select::residue_selector::SecondaryStructureSelector );
		in_selector->set_selected_ss("L");
		in_selector->set_use_dssp(true);
		protocols::pose_sewing::movers::AddFlankingVirtualResiduesMoverOP virt_mover(new protocols::pose_sewing::movers::AddFlankingVirtualResiduesMover);
		virt_mover->set_N_term_length(50);
		virt_mover->set_C_term_length(50);
		virt_mover->set_chain_to_modify(1);
		virt_mover->set_vital_selector(in_selector);
		virt_mover->add_flanking_virtual_residues(test_3_bundle);

		protocols::pose_sewing::movers::SewAnythingAddMoverOP add_mover = protocols::pose_sewing::movers::SewAnythingAddMoverOP(new protocols::pose_sewing::movers::SewAnythingAddMover());
		add_mover->set_segment_file_paths(paths);
		add_mover->set_permissible_segment_ends("H");
		add_mover->apply(test_3_bundle);
		core::select::residue_selector::VirtualResidueSelectorOP v_selector( new core::select::residue_selector::VirtualResidueSelector );
		auto returned_subset = v_selector->apply(test_3_bundle);
		core::Size real_res = 0;
		for ( core::Size cur_res = 1; cur_res <= test_3_bundle.size(); ++cur_res ) {
			if ( !returned_subset[cur_res] ) {
				++real_res;
			}
		}
		TS_ASSERT(real_res >= orig_size);
	}

	void tearDown() {
		utility::file::file_delete("test_segments/segment.txt" );
		utility::file::file_delete("test_segments/_1.pdb" );
		utility::file::file_delete("test_segments/_2.pdb" );

	}






};
