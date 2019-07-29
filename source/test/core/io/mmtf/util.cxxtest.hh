// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/mmtf/util.cxxtest.hh
/// @brief  test suite for basic mmtf helper fns
/// @author Danny Farrell

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <src/core/io/mmtf/util.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

#include <mmtf.hpp>
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>

static basic::Tracer TR("core.io.mmtf.util.cxxtest");

class mmtf_util : public CxxTest::TestSuite
{
private:
	core::pose::PoseOP pdb_pose;
	utility::vector1< std::string > mmtf_test_files;

public:
	mmtf_util() {}

	void setUp() {}

	void tearDown() {}

	void test_sd_index() {
		using namespace core::io::mmtf;
		sd_index sd_i;
		std::stringstream ss;
		ss << "mmtf::sd_index mi: " << -1 << " ci: " << -1 << " gi: " << -1 << " gai: " << -1;
		// tests == operator and string operator in one :o
		TS_ASSERT_EQUALS(sd_i.to_string(), ss.str());
	}

	void test_make_atom_num_to_sd_map() {
		using namespace core::io::mmtf;
		std::string const mmtf_filename("core/io/mmtf/1V4F.mmtf");
		::mmtf::StructureData sd;
		::mmtf::decodeFromFile(sd, mmtf_filename);
		std::map<core::Size, sd_index> const atom_num_to_sd_map(
			core::io::mmtf::make_atom_num_to_sd_map(sd));

		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(0), sd_index(0, 0, 0, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(1), sd_index(0, 0, 0, 1));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(2), sd_index(0, 0, 0, 2));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(3), sd_index(0, 0, 0, 3));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(4), sd_index(0, 0, 1, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(5), sd_index(0, 0, 1, 1));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(6), sd_index(0, 0, 1, 2));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(7), sd_index(0, 0, 1, 3));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(8), sd_index(0, 0, 1, 4));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(9), sd_index(0, 0, 1, 5));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(10), sd_index(0, 0, 1, 6));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(11), sd_index(0, 0, 2, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(12), sd_index(0, 0, 2, 1));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(13), sd_index(0, 0, 2, 2));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(14), sd_index(0, 0, 2, 3));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(15), sd_index(0, 0, 2, 4));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(16), sd_index(0, 0, 2, 5));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(17), sd_index(0, 0, 2, 6));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(18), sd_index(0, 0, 2, 7));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(19), sd_index(0, 0, 3, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(20), sd_index(0, 0, 3, 1));
		// skip a bunch
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(168), sd_index(0, 5, 56, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(169), sd_index(0, 5, 57, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(170), sd_index(0, 5, 58, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(171), sd_index(0, 5, 59, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(172), sd_index(0, 5, 60, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(173), sd_index(0, 5, 61, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(174), sd_index(0, 5, 62, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(175), sd_index(0, 5, 63, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(176), sd_index(0, 5, 64, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(177), sd_index(0, 5, 65, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(178), sd_index(0, 5, 66, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(179), sd_index(0, 5, 67, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(180), sd_index(0, 5, 68, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.at(181), sd_index(0, 5, 69, 0));
		TS_ASSERT_EQUALS(atom_num_to_sd_map.size(), (core::Size)182);
	}
};

