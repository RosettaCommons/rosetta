// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/contact_map/ContactMap.cxxtest.hh
/// @brief  test for ContactMap mover
/// @author Rocco Moretti (rmorettiase@gmail.com), Joerg Schaarschmidt (joerg.schaarschmidt@medizin.uni-leipzig.de)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/contact_map/ContactMap.hh>

// Utility Headers

static basic::Tracer TR("protocols.contact_map.ContactMap.cxxtest.hh");

// --------------- Test Class --------------- //

class ContactMapTests : public CxxTest::TestSuite {

private:
public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_contact_partner() {

		//For CA and CB the atom name isn't added, for other atoms it is.

		protocols::contact_map::ContactPartner cp(123,"HIS","CB");

		TS_ASSERT_EQUALS( cp.seqpos(), 123 );
		TS_ASSERT_EQUALS( cp.resname(), "HIS" );
		TS_ASSERT_EQUALS( cp.atomname(), "CB" );
		TS_ASSERT_EQUALS( cp.string_rep(), "HIS123" );

		protocols::contact_map::ContactPartner cp1(123,"HIS","CA");
		TS_ASSERT_EQUALS( cp1.string_rep(), "HIS123" );

		protocols::contact_map::ContactPartner cp2(123,"HIS","ND1");
		TS_ASSERT_EQUALS( cp2.string_rep(), "HIS123-ND1" );

	}

	void test_contact() {
		using namespace protocols::contact_map;

		ContactPartner cp1(123,"HIS","CB");
		ContactPartner cp2(234,"ASP","CB");

		Contact contact(cp1, cp2);

		// No counts means zero.
		TS_ASSERT_EQUALS( contact.string_rep(), "0");
		TS_ASSERT_EQUALS( contact.string_rep(4), "0.000");
		TS_ASSERT_EQUALS( contact.long_string_rep(), "HIS123\tASP234\t0");
		TS_ASSERT_EQUALS( contact.long_string_rep(4), "HIS123\tASP234\t0.000");

		// add_distance adds a contact count.
		contact.add_distance();
		TS_ASSERT_EQUALS( contact.string_rep(), "1");
		contact.add_distance();
		TS_ASSERT_EQUALS( contact.string_rep(), "2");
		contact.add_distance();
		TS_ASSERT_EQUALS( contact.string_rep(), "3");
		TS_ASSERT_EQUALS( contact.string_rep(4), "0.750");
		TS_ASSERT_EQUALS( contact.long_string_rep(), "HIS123\tASP234\t3");
		TS_ASSERT_EQUALS( contact.long_string_rep(4), "HIS123\tASP234\t0.750");

		contact.reset_count();
		TS_ASSERT_EQUALS( contact.string_rep(), "0");
		contact.add_distance();
		contact.add_distance();
		TS_ASSERT_EQUALS( contact.string_rep(), "2");

		// A numeric distance in angstroms overrides the contact count and isn't divided by nposes
		contact.add_distance(1.3454);
		TS_ASSERT_EQUALS( contact.string_rep(), "1.345");
		TS_ASSERT_EQUALS( contact.string_rep(4), "1.345");

		contact.add_distance(1.2);
		TS_ASSERT_EQUALS( contact.string_rep(), "1.200");
		TS_ASSERT_EQUALS( contact.long_string_rep(), "HIS123\tASP234\t1.200");
		TS_ASSERT_EQUALS( contact.long_string_rep(4), "HIS123\tASP234\t1.200");
	}


	void test_fill_contacts() {
		using namespace protocols::contact_map;

		core::pose::Pose pose(create_test_in_pdb_pose());

		ContactMap contact_map;
		contact_map.fill_contacts(4, 10, pose);
		//contact_map.write_to_stream( TR );

		Contact const & c1( contact_map.get_contact(1,5) );
		TS_ASSERT_EQUALS( c1.long_string_rep(), "THR4\tILE8\t0");
		Contact const & c2( contact_map.get_contact(5,1) );
		TS_ASSERT_EQUALS( c2.long_string_rep(), "THR4\tILE8\t0");
		Contact const & c3( contact_map.get_contact(3,3) );
		TS_ASSERT_EQUALS( c3.long_string_rep(), "ASP10\tASP10\t0");
	}

	void test_fill_contacts_ligand() {
		using namespace protocols::contact_map;

		core::pose::Pose pose(create_test_in_pdb_pose());

		ContactMap contact_map;
		contact_map.fill_contacts(4, 10, 29, pose);
		//contact_map.write_to_stream( TR );

		// Rows are residues in region1
		// Columns are atoms in ligand.
		Contact const & c1( contact_map.get_contact(1,1) );
		TS_ASSERT_EQUALS( c1.long_string_rep(), "THR4\tSER29- N  \t0");
		Contact const & c2( contact_map.get_contact(2,5) );
		// "space-CA-space" will get printed (CA with no spaces won't)
		TS_ASSERT_EQUALS( c2.long_string_rep(), "ILE5\tSER29- CB \t0");
		Contact const & c3( contact_map.get_contact(5,2) );
		// "space-CA-space" will get printed (CA with no spaces won't)
		TS_ASSERT_EQUALS( c3.long_string_rep(), "ILE8\tSER29- CA \t0");

		// Test if ordering matters.
		ContactMap contact_map2;
		contact_map2.fill_contacts(30, 40, 10, pose);
		Contact const & c4( contact_map2.get_contact(5,2) );
		TS_ASSERT_EQUALS( c4.long_string_rep(), "TRP34\tASP10- CA \t0");

		// Test internal ligand.
		ContactMap contact_map3;
		contact_map3.fill_contacts(30, 40, 35, pose);
		Contact const & c5( contact_map3.get_contact(3,2) );
		TS_ASSERT_EQUALS( c5.long_string_rep(), "SER32\tHIS35- CA \t0");
	}

	void test_fill_two_region() {
		using namespace protocols::contact_map;

		core::pose::Pose pose(create_test_in_pdb_pose());

		ContactMap contact_map;
		contact_map.fill_contacts(4, 10, 9, 16, pose); // intentionally overlapping
		//contact_map.write_to_stream( TR );

		Contact const & c1( contact_map.get_contact(1,1) );
		TS_ASSERT_EQUALS( c1.long_string_rep(), "THR4\tLEU9\t0");
		Contact const & c2( contact_map.get_contact(7,2) );
		TS_ASSERT_EQUALS( c2.long_string_rep(), "ASP10\tASP10\t0");
		Contact const & c3( contact_map.get_contact(2,7) );
		TS_ASSERT_EQUALS( c3.long_string_rep(), "ILE5\tASN15\t0");

		ContactMap contact_map2;
		contact_map2.fill_contacts(9, 16, 4, 10, pose); // intentionally overlapping
		//contact_map2.write_to_stream( TR );

		Contact const & c4( contact_map2.get_contact(1,1) );
		TS_ASSERT_EQUALS( c4.long_string_rep(), "LEU9\tTHR4\t0");
		Contact const & c5( contact_map2.get_contact(7,2) );
		TS_ASSERT_EQUALS( c5.long_string_rep(), "ASN15\tILE5\t0");
		Contact const & c6( contact_map2.get_contact(2,7) );
		TS_ASSERT_EQUALS( c6.long_string_rep(), "ASP10\tASP10\t0");
	}

	void test_parse_region() {
		using namespace protocols::contact_map;

		core::pose::Pose pose(create_trpcage_ideal_pose());
		ContactMap contact_map;

		core::Size begin, end;
		contact_map.parse_region_string("2-10",begin,end,pose);
		TS_ASSERT_EQUALS( begin, 2);
		TS_ASSERT_EQUALS( end, 10);

		contact_map.parse_region_string("A",begin,end,pose);
		TS_ASSERT_EQUALS( begin, 1);
		TS_ASSERT_EQUALS( end, 20);

		contact_map.parse_region_string("5",begin,end,pose);
		TS_ASSERT_EQUALS( begin, 5);
		TS_ASSERT_EQUALS( end, 5);

		contact_map.parse_region_string("10-2",begin,end,pose);
		TS_ASSERT_EQUALS( begin, 2);
		TS_ASSERT_EQUALS( end, 10);
	}

	void test_parse_tag() {
		using namespace protocols::contact_map;
		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;

		core::pose::Pose pose(create_trpcage_ideal_pose());

		ContactMap contact_map;

		TagCOP tag = tagptr_from_string("<ContactMap name=test />\n");
		contact_map.parse_my_tag( tag, data, filters, movers, pose );

		Contact const & c1( contact_map.get_contact(20,19) );
		TS_ASSERT_EQUALS( c1.long_string_rep(), "PRO19\tSER20\t0");

		ContactMap contact_map2;
		tag = tagptr_from_string("<ContactMap name=test region1=3-6/>\n");
		contact_map2.parse_my_tag( tag, data, filters, movers, pose );

		Contact const & c2( contact_map2.get_contact(2,3) );
		TS_ASSERT_EQUALS( c2.long_string_rep(), "ILE4\tGLN5\t0");

		ContactMap contact_map3;
		tag = tagptr_from_string("<ContactMap name=test region1=3-6 region2=10/>\n");
		contact_map3.parse_my_tag( tag, data, filters, movers, pose );

		Contact const & c3( contact_map3.get_contact(2,1) );
		TS_ASSERT_EQUALS( c3.long_string_rep(), "ILE4\tGLY10\t0");

		ContactMap contact_map4;
		tag = tagptr_from_string("<ContactMap name=test region1=3-6 ligand=10/>\n");
		contact_map4.parse_my_tag( tag, data, filters, movers, pose );

		Contact const & c4( contact_map4.get_contact(2,2) );
		TS_ASSERT_EQUALS( c4.long_string_rep(), "ILE4\tGLY10- CA \t0");

		ContactMap contact_map5;
		tag = tagptr_from_string("<ContactMap name=test region1=A />\n");
		contact_map5.parse_my_tag( tag, data, filters, movers, pose );

		Contact const & c5( contact_map5.get_contact(5,3) );
		TS_ASSERT_EQUALS( c5.long_string_rep(), "TYR3\tGLN5\t0");
	}

	void test_apply() {
		using namespace protocols::contact_map;
		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;

		core::pose::Pose pose(create_test_in_pdb_pose());

		ContactMap contact_map;

		TagCOP tag = tagptr_from_string("<ContactMap name=test region1=4-10 models_per_file=0 />\n"); //don't output textfile
		contact_map.parse_my_tag( tag, data, filters, movers, pose );

		contact_map.apply(pose);
		TS_ASSERT_EQUALS( contact_map.get_contact(1,7).long_string_rep(), "THR4\tASP10\t0"); // 11.11
		TS_ASSERT_EQUALS( contact_map.get_contact(2,5).long_string_rep(), "ILE5\tILE8\t1"); // 5.75

		contact_map.apply(pose);
		TS_ASSERT_EQUALS( contact_map.get_contact(1,7).long_string_rep(), "THR4\tASP10\t0"); // 11.11
		TS_ASSERT_EQUALS( contact_map.get_contact(2,5).long_string_rep(), "ILE5\tILE8\t2"); // 5.75

		// Test distance cutoff
		ContactMap contact_map2;

		tag = tagptr_from_string("<ContactMap name=test region1=4-10 models_per_file=0 distance_cutoff=15 />\n");
		contact_map2.parse_my_tag( tag, data, filters, movers, pose );

		contact_map2.apply(pose);
		contact_map2.apply(pose);
		TS_ASSERT_EQUALS( contact_map2.get_contact(1,7).long_string_rep(), "THR4\tASP10\t2"); // 11.11
		TS_ASSERT_EQUALS( contact_map2.get_contact(2,5).long_string_rep(), "ILE5\tILE8\t2"); // 5.75

		// Distance matrix form
		ContactMap contact_map3;

		tag = tagptr_from_string("<ContactMap name=test region1=4-10 models_per_file=0 distance_matrix=true/>\n"); //don't output textfile
		contact_map3.parse_my_tag( tag, data, filters, movers, pose );

		contact_map3.apply(pose);
		TS_ASSERT_EQUALS( contact_map3.get_contact(1,7).long_string_rep(), "THR4\tASP10\t11.108"); // 11.108
		TS_ASSERT_EQUALS( contact_map3.get_contact(2,5).long_string_rep(), "ILE5\tILE8\t5.752"); // 5.752
	}

};
