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
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>

#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ResidueSubGraphTypes.hh>

// Project Headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/pose/Pose.hh>
// Platform Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <ostream>

#include <core/import_pose/import_pose.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <core/chemical/AtomTypeSet.fwd.hh>

#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>


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

};
