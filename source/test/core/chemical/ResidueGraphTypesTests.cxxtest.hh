// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueGraphTypeTests.cxxtest.hh
/// @brief unit tests for ResidueType
/// @author Steven Combs


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/MutableResidueType.hh>

#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ResidueSubGraphTypes.hh>

// Project Headers
#include <core/pose/Pose.hh>
// Platform Headers

// C++ Headers
#include <string>

#include <core/chemical/Atom.hh>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/filtered_graph.hpp>

#include <test/util/pose_funcs.hh>


class ResidueGraphTypesTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_heavy_atom_filter() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(1));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::HeavyAtomGraph HAgraph( core::chemical::make_heavy_atom_graph( rsd ) );

		core::chemical::HeavyAtomVIterPair vp = boost::vertices(HAgraph);
		core::chemical::VD vd = *vp.first; //get the vertex

		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), " N  ");
		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), " CA ");
	}

	void test_acceptor_atom_filter() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(1));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::AcceptorAtomGraph HAgraph( core::chemical::make_acceptor_atom_graph( rsd ) );

		core::chemical::AcceptorAtomVIterPair vp = boost::vertices(HAgraph);
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), " O  ");
		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), " OD1");


	}

	void test_heavy_atom_with_polar_hydrogen_filter(){
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(1));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::HeavyAtomWithPolarHydrogensGraph HAgraph( core::chemical::make_heavy_atom_with_polar_hydrogens_graph( rsd ) );
		core::chemical::HeavyAtomWithPolarHydrogensVIterPair vp = boost::vertices(HAgraph);
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), " N  ");
		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), " ND2");


	}

	void test_heavy_atom_with_hydrogens_filter(){
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(1));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::HeavyAtomWithHydrogensGraph HAgraph( core::chemical::make_heavy_atom_with_hydrogens_graph( rsd ) );
		core::chemical::HeavyAtomWithHydrogensVIterPair vp = boost::vertices(HAgraph);
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), " N  ");
		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), " CA ");
	}

	void test_hydrogen_atom_filter() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(1));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::HydrogenAtomGraph HAgraph( core::chemical::make_hydrogen_atom_graph( rsd ) );

		//now lets iterate through the graph
		std::set< std::string > atom_names; // We don't care about ordering
		core::chemical::HydrogenAtomVIterPair vp = boost::vertices(HAgraph);
		for ( ; vp.first != vp.second; ++vp.first ) {
			atom_names.insert( HAgraph[*vp.first].name() );
		}
		std::set< std::string > reference{ "1H  ", "2H  ", "3H  ", " HA ", "1HB ", "2HB ", "1HD2", "2HD2" };
		TS_ASSERT_EQUALS( atom_names, reference );
	}


	void test_polar_hydrogen_atom_filter() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(1));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::PolarHydrogenGraph HAgraph( core::chemical::make_polar_hydrogen_graph( rsd ) );

		//now lets iterate through the graph
		std::set< std::string > atom_names; // We don't care about ordering
		core::chemical::PolarHydrogenVIterPair vp = boost::vertices(HAgraph);
		for ( ; vp.first != vp.second; ++vp.first ) {
			atom_names.insert( HAgraph[*vp.first].name() );
		}
		std::set< std::string > reference{ "1H  ", "2H  ", "3H  ", "1HD2", "2HD2" };
		TS_ASSERT_EQUALS( atom_names, reference );
	}

	void test_aromatic_hydrogen_filter(){
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(3));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::AromaticAtomGraph HAgraph( core::chemical::make_aromatic_atom_graph( rsd ) );

		core::chemical::AromaticAtomVIterPair vp = boost::vertices(HAgraph);
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), /*" CG "*/" CD1");

		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), /*" CD1"*/" CD2");

	}

	void test_apolar_hydrogen_atom_filter() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType const & rsd_origin(pose.residue_type(1));
		core::chemical::MutableResidueType rsd( rsd_origin );
		core::chemical::APolarHydrogenGraph HAgraph( core::chemical::make_apolar_hydrogen_graph( rsd ) );

		//now lets iterate through the heavy atom graph
		core::chemical::APolarHydrogenVIterPair vp = boost::vertices(HAgraph);
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), " HA ");
		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), "1HB ");


	}

	void test_vertex_index(){
		core::chemical::ResidueGraph graph;
		utility::vector0<core::chemical::VD> storage;
		for ( core::Size i=0; i<= 9; ++i ) {
			storage.push_back( graph.add_vertex());
		}

		for ( core::Size i=0; i<= 9; ++i ) {
			TS_ASSERT_EQUALS(i, boost::get(boost::vertex_index, graph, storage[i]));
		}

		graph.remove_vertex(storage[5]);
		graph.remove_vertex(storage[6]);

		core::chemical::regenerate_graph_vertex_index(graph);
		core::chemical::VIter itr, end;

		core::Size count(0);
		for ( boost::tie(itr,end) = boost::vertices(graph); itr !=end; ++itr ) {
			TS_ASSERT_EQUALS(count, boost::get(boost::vertex_index, graph, *itr));
			++count;
			//boost::put(boost::vertex_index, graph, *itr, index);
		}

		core::chemical::VD vd = graph.add_vertex();
		core::chemical::regenerate_graph_vertex_index(graph);

		TS_ASSERT_EQUALS(boost::get(boost::vertex_index, graph, vd), 8);


	}

};
