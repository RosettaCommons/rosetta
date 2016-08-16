// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/AtomTypeSet.cxxtest.hh

/// @brief  test suite for reading and writing atom type sets to and from a database
/// @author Matthew O'Meara (mattjomeara@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeDatabaseIO.hh>
#include <basic/database/open.hh>
#include <basic/database/sql_utils.hh>
#include <core/types.hh>

// Project headers
#include <basic/Tracer.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <map>

static basic::Tracer tr("core.chemical.AtomTypeDatabaseIOTests.cxxtest");

class AtomTypeDatabaseIOTests : public CxxTest::TestSuite {
public:

	AtomTypeDatabaseIOTests()
	{
		core_init();

		std::string database_filename("atom_type_set_test_database.db3");
		utility::file::file_delete(database_filename);
		db_session_ = basic::database::get_db_session(database_filename);
		atom_type_dbio_.initialize(db_session_);
	}

	virtual ~AtomTypeDatabaseIOTests() {}

	static AtomTypeDatabaseIOTests *createSuite() {
		return new AtomTypeDatabaseIOTests();
	}

	static void destroySuite( AtomTypeDatabaseIOTests *suite ) {
		delete suite;
	}


	// --------------- Test Cases --------------- //

	void test_atom_type_set_name() {
		using basic::database::full_name;
		using core::Size;
		using core::chemical::AtomType;
		using core::chemical::AtomTypeSet;
		using core::chemical::AtomTypeSetOP;
		using std::string;
		using std::map;
		using utility::vector1;

		vector1<string> ats_directories;
		ats_directories.push_back(
			full_name("chemical/atom_type_sets/centroid"));
		ats_directories.push_back(
			full_name("chemical/atom_type_sets/centroid_min"));
		ats_directories.push_back(
			full_name("chemical/atom_type_sets/coarse_rna"));
		ats_directories.push_back(
			full_name("chemical/atom_type_sets/coarse_two_bead"));
		ats_directories.push_back(
			full_name("chemical/atom_type_sets/fa_standard"));

		// create atom type sets from the rosetta_database
		vector1<AtomTypeSetOP> atss;  // "Atom Type SetS"
		foreach ( string ats_directory, ats_directories ) {
			atss.push_back(AtomTypeSetOP( new AtomTypeSet(ats_directory) ));
		}

		// save atom type set to ats_test_database.db3
		foreach ( AtomTypeSetOP ats, atss ) {
			atom_type_dbio_.write_atom_type_set_to_database(*ats, db_session_);
		}

		// create atom type sets from ats_test_database.db3
		vector1<AtomTypeSetOP> atss_from_db;
		foreach ( AtomTypeSetOP ats, atss ) {
			string const & ats_name(ats->name());
			atss_from_db.push_back(AtomTypeSetOP( new AtomTypeSet(ats_name, db_session_) ));
		}

		// test that the atom type set has been preserved after reading
		// and writing it to atss_test_database.db3
		for ( Size i=1; i <= atss.size(); ++i ) {
			AtomTypeSet const & ats_orig(*atss[i]);
			AtomTypeSet const & ats_new(*atss_from_db[i]);

			TS_ASSERT_EQUALS(ats_orig.name(), ats_new.name());
			TS_ASSERT_EQUALS(ats_orig.n_atomtypes(), ats_new.n_atomtypes());
			TS_ASSERT_EQUALS(ats_orig.directory(), ats_new.directory());

			TS_ASSERT_EQUALS(
				ats_orig.extra_parameter_indices().size(),
				ats_new.extra_parameter_indices().size());

			// check the atom types contain the same information
			for ( Size i=1; i <= ats_orig.n_atomtypes(); ++i ) {
				locate_matching_atom_type(ats_orig, i, ats_new);
			}
		}
	}


	// Assert that the there is an atom type in ats_new that matches
	// atom type i in ats_orig
	void
	locate_matching_atom_type(
		core::chemical::AtomTypeSet const & ats_orig,
		core::Size const i,
		core::chemical::AtomTypeSet const & ats_new
	) const {
		core::chemical::AtomType const & a(ats_orig[i]);
		bool found_match(false);
		for ( core::Size j=1; j <= ats_new.n_atomtypes(); ++j ) {
			core::chemical::AtomType const & b(ats_new[j]);
			if ( !a.name().compare(b.name()) ) {
				found_match = true;
				TS_ASSERT_EQUALS(a.element(), b.element());
				TS_ASSERT_EQUALS(a.lj_radius(), b.lj_radius());
				TS_ASSERT_EQUALS(a.lj_wdepth(), b.lj_wdepth());
				TS_ASSERT_EQUALS(a.lk_lambda(), b.lk_lambda());
				TS_ASSERT_EQUALS(a.lk_volume(), b.lk_volume());
				TS_ASSERT_EQUALS(a.lk_dgfree(), b.lk_dgfree());
				TS_ASSERT_EQUALS(a.is_acceptor(), b.is_acceptor());
				TS_ASSERT_EQUALS(a.is_donor(), b.is_donor());
				TS_ASSERT_EQUALS(a.is_polar_hydrogen(), b.is_polar_hydrogen());
				TS_ASSERT_EQUALS(a.is_h2o(), b.is_h2o());
				TS_ASSERT_EQUALS(a.is_aromatic(), b.is_aromatic());
				TS_ASSERT_EQUALS(a.atom_has_orbital(), b.atom_has_orbital());
				TS_ASSERT_EQUALS(a.is_virtual(), b.is_virtual());
				TS_ASSERT_EQUALS(a.hybridization(), b.hybridization());

				for ( std::map<std::string, int>::const_iterator
						extra_parameter_index_iter =
						ats_orig.extra_parameter_indices().begin(),
						extra_parameter_index_iter_end =
						ats_orig.extra_parameter_indices().end();
						extra_parameter_index_iter != extra_parameter_index_iter_end;
						++extra_parameter_index_iter ) {
					std::string const & param(extra_parameter_index_iter->first);
					int a_index(extra_parameter_index_iter->second);
					int b_index(ats_new.extra_parameter_index(param));
					TS_ASSERT_EQUALS(a.extra_parameter(a_index), b.extra_parameter(b_index));
				}
			}
		}
		TS_ASSERT(found_match);
	}

private:
	core::chemical::AtomTypeDatabaseIO atom_type_dbio_;
	utility::sql_database::sessionOP db_session_;
};
