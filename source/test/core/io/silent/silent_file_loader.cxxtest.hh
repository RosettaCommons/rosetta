// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/silent_file_loader.cxxtest.hh
/// @brief test suite for core::io::silent::SilentFileLoader
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileLoader.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <basic/resource_manager/types.hh>

// Platform headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <algorithm>


static basic::Tracer TR("core.io.silent.silent_file_loader");


class SilentFileLoaderTests : public CxxTest::TestSuite {

public:

	void setUp() {
				core_init_with_additional_options( "-mute core.io.pdb -mute core.conformation" );
	}

	void test_all(){
		do_test_basic(1e-2, 1e-1);
	}

	void do_test_basic(
		double const rms_threshold,
		double const score_threshold
	) {

		core::io::silent::SilentFileLoader sfl;
		core::io::silent::SilentFileOptions sfo;
		try {
			sfo.set_silent_struct_type("BLA");
			TS_ASSERT( false );
		} catch (utility::excn::EXCN_BadInput & e){
			TR << e.msg();
		}

		sfo.set_silent_struct_type("binary");

		// SETUP THE STARTING SILENT FILE
		core::pose::Pose ref_pose, restored_pose;
		core::chemical::ResidueTypeSetCAP rsd =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::import_pose::pose_from_file( ref_pose, *rsd, std::string("core/io/bin_silentfile_test.pdb"), core::import_pose::PDB_file);
		std::string const silent_outfile( "core/io/bin_silentfile_test.out" ); // read file w/ non-ideal geometry
		utility::file::file_delete( silent_outfile );
		core::io::silent::SilentFileData sfd;
		core::io::silent::BinaryProteinSilentStruct pss( ref_pose, "tag" );
		sfd.write_silent_struct( pss, silent_outfile );


		// READ IN SILENT FILE
		utility::io::izstream data( silent_outfile.c_str() );
		basic::resource_manager::ResourceOP resource(
			sfl.create_resource(sfo, "locator_id", data));
		restored_pose = static_cast<core::pose::Pose const &>(*resource);

		// TEST RMS DIFFERENCE
		core::Real rms_to_restored = core::scoring::CA_rmsd( ref_pose, restored_pose );
		TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test score13 difference
		core::scoring::ScoreFunctionOP scorefxn =
					core::scoring::ScoreFunctionFactory::create_score_function( "score13_env_hb" );
		core::Real score_ref = (*scorefxn)(ref_pose);
		core::Real score_restored = (*scorefxn)(restored_pose);
		core::Real score_del = std::fabs( score_restored - score_ref );
		TR << "Score difference: " << score_del << std::endl;
		if ( score_del > score_threshold ) {
			restored_pose.dump_pdb( "restored_pose_wrong_score.pdb" );
			ref_pose.dump_pdb( "ref_pose_wrong_score.pdb" );
		}
		TS_ASSERT( score_del < score_threshold );

	}
};
