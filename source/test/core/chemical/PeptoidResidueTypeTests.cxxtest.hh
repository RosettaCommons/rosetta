// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/chemical/PeptoidResidueTypeTests.cxxtest.hh
/// @brief  Unit tests to confirm that we can build reasonable peptoid geometry.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/id/NamedAtomID.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>

static basic::Tracer TR("PeptoidResidueTypeTests");

#define OMEGA_EPSILON_CHIRAL 3
#define OMEGA_EPSILON 2

class PeptoidResidueTypeTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	core::pose::PoseOP build_peptoid_pose() {
		std::stringstream sequence_accumulator;
		sequence_accumulator << "G"; //Pad with glycine.

		for ( core::Size i(1), imax(999); i<=imax; ++i ) { //Try to make one peptoid of each possible type.
			std::string const candidate_name( std::to_string(i) );
			core::chemical::ResidueTypeFinder finder( *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) );
			core::chemical::ResidueTypeCOP restype( finder.residue_base_name( candidate_name ).get_representative_type(false) );
			if ( restype == nullptr ) continue;
			sequence_accumulator << "X[" << candidate_name << "]";
			TS_ASSERT( restype->is_peptoid() );
			if ( i >= 600 && i < 700 ) {
				TS_ASSERT( restype->is_r_peptoid() || restype->is_s_peptoid() );
				if ( i % 2 == 0 ) {
					TS_ASSERT( restype->is_r_peptoid() );
				} else {
					TS_ASSERT( restype->is_s_peptoid() );
				}
			} else {
				TS_ASSERT( !( restype->is_r_peptoid() || restype->is_s_peptoid() ) );
			}
		}

		if ( sequence_accumulator.str() == "G" ) return nullptr;
		sequence_accumulator << "G"; //Pad with glycine.

		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
		core::pose::make_pose_from_sequence( *pose, sequence_accumulator.str(), core::chemical::FA_STANDARD, true, false );

		for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
			pose->set_phi( ir, 170 );
			pose->set_psi( ir, 170 );
			pose->set_omega( ir, 180 );

			if ( ir == 1 || ir == irmax ) continue;

			// Check that we have the right atom names:
			core::chemical::ResidueType const & restype( pose->residue_type(ir) );
			TS_ASSERT( restype.has( "N" ) );
			TS_ASSERT( restype.has( "CA" ) );
			TS_ASSERT( restype.has( "C" ) );
			TS_ASSERT( restype.has( "O" ) );
			TS_ASSERT( restype.has( "CA1" ) );
		}
		pose->update_residue_neighbors();

		return pose;
	}


	void test_lower_connect_dependence(){
		core::pose::PoseOP pose( build_peptoid_pose() );
		TS_ASSERT( pose != nullptr );
		TS_ASSERT( pose->total_residue() > 2 ); //If length == 2, then we just have the terminal glycines.

		core::Size badcount(0);

		for ( core::Size i(2), imax(pose->total_residue()); i<imax; ++i ) {
			core::Real const epsilon( ( pose->residue_type(i).is_s_peptoid() || pose->residue_type(i).is_r_peptoid() ) ? OMEGA_EPSILON_CHIRAL : OMEGA_EPSILON );

			core::Real angle1( 0.0 ), angle2( 0.0 );
			numeric::dihedral_degrees_double( pose->xyz( core::id::NamedAtomID( "O", i-1 ) ), pose->xyz( core::id::NamedAtomID( "C", i-1 ) ), pose->xyz( core::id::NamedAtomID( "N", i ) ), pose->xyz( core::id::NamedAtomID( "CA1", i ) ), angle1 );
			angle1 = numeric::nonnegative_principal_angle_degrees( angle1 );

			pose->set_omega( i-1, 77.3 );
			numeric::dihedral_degrees_double( pose->xyz( core::id::NamedAtomID( "O", i-1 ) ), pose->xyz( core::id::NamedAtomID( "C", i-1 ) ), pose->xyz( core::id::NamedAtomID( "N", i ) ), pose->xyz( core::id::NamedAtomID( "CA1", i ) ), angle2 );
			angle2 = numeric::nonnegative_principal_angle_degrees( angle2 );

			TR << "Peptoid " << pose->residue_type(i).base_name() << ": angle1=" << angle1 << ", angle2=" << angle2;
			if ( std::abs( angle1 - 180 ) > epsilon || std::abs( angle2 - 77.3 ) > epsilon ) {
				TR << " <--------------- BAD PEPTOID GEOMETRY!!!";
				++badcount;
			}
			TR << std::endl;
			TS_ASSERT_DELTA( angle1, 180, epsilon );
			TS_ASSERT_DELTA( angle2, 77.3, epsilon );
		}
		if ( badcount != 0 ) {
			TR << badcount << " peptoid residues produce bad geometry because sidechains don't move properly when omega changes." << std::endl;
		} else {
			TR << "All peptoid residues have sidechains that move properly when omega changes." << std::endl;
		}
		//pose->dump_pdb( "PEPTOID_TEMP.pdb" ); //DELETE ME
	}

};
