// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/AtomICoor.cxxtest.hh
/// @brief unit tests for AtomICoor
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <core/init_util.hh>

// Unit Headers
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Platform Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <ostream>

//Auto Headers

static basic::Tracer TR("core.chemical.AtomICoor.cxxtest");

class AtomICoorTests : public CxxTest::TestSuite {

	core::pose::PoseCOP pose_;

public:

	void setUp() {
		core_init_with_additional_options("-extra_res_fa core/chemical/ASX.params core/chemical/LYX.params"); // Remove if database versions get permanently enabled
		pose_ = core::import_pose::pose_from_file("core/chemical/conn.pdb");
	}

	void tearDown() {}

	void test_atom_icoord_build() {
		core::conformation::Residue const & rsd( pose_->residue(2) );
		TS_ASSERT_EQUALS( rsd.name(), "ASX" );

		core::chemical::AtomICoor const & ha( rsd.icoor( rsd.atom_index("HA") ) );

		core::Vector pos( ha.build( rsd, pose_->conformation() ) );

		TS_ASSERT_DELTA( pos.x(), rsd.xyz("HA").x(), 0.001 );
		TS_ASSERT_DELTA( pos.y(), rsd.xyz("HA").y(), 0.001 );
		TS_ASSERT_DELTA( pos.z(), rsd.xyz("HA").z(), 0.001 );
		TR << "test_atom_icoord_build Done!" << std::endl;
	}

	void test_id_poly() {
		core::conformation::Residue const & rsd( pose_->residue(2) );

		// Internal
		core::chemical::AtomICoor const & ha( rsd.icoor( rsd.atom_index("HA") ) );

		core::id::AtomID ha_stub1_id( ha.stub_atom1().atom_id( rsd, pose_->conformation() ) );
		TS_ASSERT_EQUALS( ha_stub1_id.rsd(), 2 );
		TS_ASSERT_EQUALS( ha_stub1_id.atomno(), rsd.atom_index("CA") );

		core::id::AtomID ha_stub3_id( ha.stub_atom3().atom_id( rsd, pose_->conformation() ) );
		TS_ASSERT_EQUALS( ha_stub3_id.rsd(), 2 );
		TS_ASSERT_EQUALS( ha_stub3_id.atomno(), rsd.atom_index("CB") );

		// Polymer Lower
		core::chemical::AtomICoor const & bbh( rsd.icoor( rsd.atom_index("H") ) );

		core::id::AtomID bbh_stub3_id( bbh.stub_atom3().atom_id( rsd, pose_->conformation() ) );
		TS_ASSERT_EQUALS( bbh_stub3_id.rsd(), 1 );
		TS_ASSERT_EQUALS( bbh_stub3_id.atomno(), pose_->residue(1).atom_index("C") );

		// Polymer Upper
		core::chemical::AtomICoor const & bbo( rsd.icoor( rsd.atom_index("O") ) );

		core::id::AtomID bbo_stub3_id( bbo.stub_atom3().atom_id( rsd, pose_->conformation() ) );
		TS_ASSERT_EQUALS( bbo_stub3_id.rsd(), 3 );
		TS_ASSERT_EQUALS( bbo_stub3_id.atomno(), pose_->residue(3).atom_index("N") );

		TR << "test_if_poly Done!" << std::endl;
	}

	void test_id_nonpoly() {
		// Test what happens with connection IDs when we delete connections.
		core::pose::Pose pose(*pose_); // make copy
		core::conformation::Residue old_res( pose.residue(2) );

		core::chemical::ResidueTypeSetCOP rt_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		// First, some sanity checks to make sure the connection information is set up correctly.

		TS_ASSERT_EQUALS( pose.residue(2).n_current_residue_connections(), 3);
		TS_ASSERT_EQUALS( pose.residue(5).n_current_residue_connections(), 3);
		TS_ASSERT_EQUALS( pose.residue(2).n_non_polymeric_residue_connections(), 1);
		TS_ASSERT_EQUALS( pose.residue(5).n_non_polymeric_residue_connections(), 1);
		// Connection 3 should be the sidechain connection
		TS_ASSERT_EQUALS( pose.residue(2).residue_connection_partner(3), 5 );
		TS_ASSERT_EQUALS( pose.residue(5).residue_connection_partner(3), 2 );
		// connect_atom is the atom on *this* residue that's connected to the parameter residue
		TS_ASSERT_EQUALS( pose.residue(2).connect_atom(pose.residue(5)), pose.residue(2).atom_index("CG") );
		TS_ASSERT_EQUALS( pose.residue(5).connect_atom(pose.residue(2)), pose.residue(5).atom_index("NZ") );

		// Now test the id for the connection.
		{
			core::chemical::AtomICoor const & sco( pose.residue(2).icoor( pose.residue(2).atom_index("OD1") ) );

			core::id::AtomID sco_stub3_id( sco.stub_atom3().atom_id( pose.residue(2), pose_->conformation() ) );
			TS_ASSERT_EQUALS( sco_stub3_id.rsd(), 5 );
			TS_ASSERT_EQUALS( pose.residue(sco_stub3_id.rsd()).atom_name( sco_stub3_id.atomno() ), " NZ " );
		}

		// Remove c-term connection.
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::UPPER_TERMINUS_VARIANT, 2 );
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::LOWER_TERMINUS_VARIANT, 3 );

		TS_ASSERT( core::pose::is_upper_terminus(pose, 2 ) );

		{
			core::chemical::AtomICoor const & sco( pose.residue(2).icoor( pose.residue(2).atom_index("OD1") ) );

			core::id::AtomID sco_stub3_id( sco.stub_atom3().atom_id( pose.residue(2), pose_->conformation() ) );
			TS_ASSERT_EQUALS( sco_stub3_id.rsd(), 5 );
			TS_ASSERT_EQUALS( pose.residue(sco_stub3_id.rsd()).atom_name( sco_stub3_id.atomno() ), " NZ " );
		}

		// Also remove the n-term connection
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::UPPER_TERMINUS_VARIANT, 1 );
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::LOWER_TERMINUS_VARIANT, 2 );

		TS_ASSERT( core::pose::is_upper_terminus(pose, 2 ) );
		TS_ASSERT( core::pose::is_lower_terminus(pose, 2 ) );

		{
			core::chemical::AtomICoor const & sco( pose.residue(2).icoor( pose.residue(2).atom_index("OD1") ) );

			core::id::AtomID sco_stub3_id( sco.stub_atom3().atom_id( pose.residue(2), pose_->conformation() ) );
			TS_ASSERT_EQUALS( sco_stub3_id.rsd(), 5 );
			TS_ASSERT_EQUALS( pose.residue(sco_stub3_id.rsd()).atom_name( sco_stub3_id.atomno() ), " NZ " );
		}

		TR << "test_id_nonpoly Done!" << std::endl;
	}

};
