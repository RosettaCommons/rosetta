// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/io/pdb/chains_62.cxxtest.hh
/// @brief   PDB chain legal chars are [A-Za-z0-9].  This tests that Rosetta knows that.
/// @author  Steven Lewis <smlewi@gmail.com>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit header

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Basic header

#include <basic/Tracer.hh>

// C++ header

static basic::Tracer TR("core.io.chains_62.cxxtest");

inline
core::pose::Pose pose_maker_62() {

	std::string const pdbstring(
		"HETATM  427 ZN    ZN a 101      69.086  11.890  98.667  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN b 101       4.071  38.062  31.605  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN c 101      61.515   0.456  40.886  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN d 101      91.445  91.670  50.968  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN e 101      16.647  70.591  33.778  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN f 101      24.806  37.393  72.569  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN g 101      17.081  37.458   8.600  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN h 101      14.436  74.413  19.489  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN i 101      78.247  37.932  76.366  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN j 101      46.036  47.337   5.197  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN k 101      32.757  27.623  51.265  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN l 101      26.094  33.680  48.024  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN m 101      32.324  78.233  97.909  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN n 101       5.753  27.577  72.453  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN o 101      46.786  95.034  71.271  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN p 101      27.607  64.235  10.943  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN q 101      94.740  70.529  80.121  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN r 101      57.596  73.329  96.508  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN s 101      22.505  51.481  28.693  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN t 101      71.655  89.577   6.724  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN u 101      65.919  21.956  18.974  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN v 101      92.464  74.266   6.619  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN w 101       1.133  93.836  55.489  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN x 101      26.209  47.353  67.372  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN y 101      80.373  70.940  77.568  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN z 101      66.518  39.232  31.299  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN A 101      77.507   7.930  47.300  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN B 101      15.816  84.927  52.088  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN C 101       9.759  67.021  74.104  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN D 101      25.623  67.339  70.821  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN E 101      77.672   8.909  59.495  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN F 101      35.843  87.611  65.519  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN G 101      71.038  53.055  92.986  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN H 101      86.081  77.328   1.863  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN I 101       1.322  98.742  42.142  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN J 101      51.759  11.421  59.208  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN K 101      36.431  45.106  31.480  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN L 101       8.001  59.100  88.468  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN M 101      36.057  65.531  54.182  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN N 101      49.203   2.511  96.125  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN O 101      24.885  80.228  35.672  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN P 101      27.809  36.885  93.313  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN Q 101      91.204  23.266  26.225  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN R 101      88.194  91.142  17.085  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN S 101       1.468  40.669  22.645  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN T 101      44.041  26.339  62.676  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN U 101      57.901  62.857   6.386  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN V 101      31.777  32.898  34.365  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN W 101      37.949  56.718  15.601  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN X 101      30.216  39.407  67.670  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN Y 101      91.769   9.261  54.108  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN Z 101      67.254  87.783  77.248  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 0 101      86.902  49.578  93.000  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 1 101      73.724  59.340  63.248  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 2 101      77.264  65.005  11.320  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 3 101      73.524  63.984  89.795  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 4 101      36.441  98.545  19.566  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 5 101      30.633  98.547  84.038  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 6 101      97.990  73.174  33.627  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 7 101      45.242  32.643  24.300  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 8 101      35.358  86.110   0.231  0.86  7.74          ZN  \n"
		"HETATM  427 ZN    ZN 9 101      28.401  63.642  66.749  0.86  7.74          ZN  \n"
	);

	return fullatom_pose_from_string(pdbstring);
}

class Chains62Tests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////

	/// @details PDB chain legal chars are [A-Za-z0-9].  Read in such a PDB, and validate that the Pose has 62 chains labeled correctly
	void test_read_62_chain_pdb()
	{

		core::pose::Pose const pose(pose_maker_62());

		core::Size const num_chains(62);
		TS_ASSERT_EQUALS(pose.num_chains(), 62);

		utility::vector1<char> chains = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

		//assert(chains.size() == num_chains);

		for ( core::Size i(1); i<=num_chains; ++i ) {
			TS_ASSERT_EQUALS(pose.pdb_info()->chain(i), chains[i]);
		}
	}
};  // class Chains62Tests
