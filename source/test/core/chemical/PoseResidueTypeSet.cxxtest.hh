// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/chemical/PoseResidueTypeSet.cxxtest.hh
/// @brief  Test ResidueTypeSetPose
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("core.chemical.PoseResidueTypeSet.cxxtest");


class PoseResidueTypeSetTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}


	void test_addition() {
		using namespace core::chemical;

		GlobalResidueTypeSetOP grs( new GlobalResidueTypeSet( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ) ) );

		ResidueType const & serine = grs->name_map( "SER" );

		ResidueTypeOP modser = serine.clone();
		modser->nbr_radius( 15.0);
		modser->name( "bigser" );

		PoseResidueTypeSetOP rs( new PoseResidueTypeSet( grs ) );

		//get some stuff from the residue type set
		core::Size n_base_res_types   = rs->base_residue_types().size();
		core::Size n_unpatchable_res_types = rs->unpatchable_residue_types().size();
		core::Size n_ser_types = ResidueTypeFinder( *rs ).name3( "SER" ).get_all_possible_residue_types().size();
		core::Size n_gln_types = ResidueTypeFinder( *rs ).name3( "GLN" ).get_all_possible_residue_types().size();
		core::Size n_ser_aa = ResidueTypeFinder( *rs ).aa( aa_ser ).get_all_possible_residue_types().size();

		//now change the residue type set
		// const_cast, for testing purposes only.
		rs->add_unpatchable_residue_type( modser );

		//now make sure everything is as should be
		TS_ASSERT( n_base_res_types == rs->base_residue_types().size());
		TS_ASSERT( n_unpatchable_res_types + 1 == rs->unpatchable_residue_types().size());
		TS_ASSERT( n_ser_types + 1 == ResidueTypeFinder( *rs ).name3( "SER" ).get_all_possible_residue_types().size() );
		TS_ASSERT( n_gln_types == ResidueTypeFinder( *rs ).name3( "GLN" ).get_all_possible_residue_types().size() );
		TS_ASSERT( n_ser_aa + 1 == ResidueTypeFinder( *rs ).aa( aa_ser ).get_all_possible_residue_types().size() );
		TS_ASSERT( rs->has_name("bigser") );
	}





};



