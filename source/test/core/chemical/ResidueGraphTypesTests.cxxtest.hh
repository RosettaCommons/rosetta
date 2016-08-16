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

#include <core/chemical/ResidueGraphTypes.hh>

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
		core::chemical::ResidueType rsd(pose.residue_type(1));
		core::chemical::HeavyAtomGraph HAgraph(rsd.heavy_atoms());

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
		core::chemical::ResidueType rsd(pose.residue_type(1));
		core::chemical::AcceptorAtomGraph HAgraph(rsd.acceptor_atoms());

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
		core::chemical::ResidueType rsd(pose.residue_type(1));
		core::chemical::HeavyAtomWithPolarHydrogensGraph HAgraph(rsd.heavy_atom_with_polar_hydrogens());
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
		core::chemical::ResidueType rsd(pose.residue_type(1));
		core::chemical::HeavyAtomWithHydrogensGraph HAgraph(rsd.heavy_atom_with_hydrogens());
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
		core::chemical::ResidueType rsd(pose.residue_type(1));
		core::chemical::HydrogenAtomGraph HAgraph(rsd.hydrogens());

		//now lets iterate through the heavy atom graph
		core::chemical::HydrogenAtomVIterPair vp = boost::vertices(HAgraph);
		//core::chemical::VIter v_iter=  //creat an iterator for access
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), "1HD2");
		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), "2HD2");

	}


	void test_polar_hydrogen_atom_filter() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType rsd(pose.residue_type(1));
		core::chemical::PolarHydrogenGraph HAgraph(rsd.polar_hydrogens());

		core::chemical::PolarHydrogenVIterPair vp = boost::vertices(HAgraph);
		//core::chemical::VIter v_iter=  //creat an iterator for access
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), "1HD2");
		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), "2HD2");

	}

	void test_aromatic_hydrogen_filter(){
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType rsd(pose.residue_type(3));
		core::chemical::AromaticAtomGraph HAgraph(rsd.aromatic_atoms());

		core::chemical::AromaticAtomVIterPair vp = boost::vertices(HAgraph);
		core::chemical::VD vd = *vp.first; //get the vertex
		core::chemical::Atom atom(HAgraph[vd]); //get the atom
		TS_ASSERT_EQUALS(atom.name(), " CG ");

		++vp.first;
		vd = *vp.first;
		atom = HAgraph[vd];
		TS_ASSERT_EQUALS(atom.name(), " CD1");


	}

	void test_apolar_hydrogen_atom_filter() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType rsd(pose.residue_type(1));
		core::chemical::APolarHydrogenGraph HAgraph(rsd.apolar_hydrogens());

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
