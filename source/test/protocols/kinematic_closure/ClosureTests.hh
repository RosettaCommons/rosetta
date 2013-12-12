// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Coordinate information for the closure tests.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_test_protocols_kinematic_closure_ClosureTests_HH
#define INCLUDED_test_protocols_kinematic_closure_ClosureTests_HH

#include <test/protocols/kinematic_closure/TestHelpers.hh>
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.hh>

using core::PointPosition;

class FiveResidueTest : public TestHelper { // {{{1
public:
	FiveResidueTest() : TestHelper("five", 5, 2, 3, 4, 2) {
		atom_xyzs[1]  = PointPosition(29.426, 38.43, 15.446);
		atom_xyzs[2]  = PointPosition(30.225, 38.643, 16.662);
		atom_xyzs[3]  = PointPosition(29.664, 39.839, 17.434);
		atom_xyzs[4]  = PointPosition(30.132, 40.069, 18.642);
		atom_xyzs[5]  = PointPosition(29.607, 41.18, 19.467);
		atom_xyzs[6]  = PointPosition(28.234, 41.641, 19.021);
		atom_xyzs[7]  = PointPosition(28.161, 42.43, 17.933);
		atom_xyzs[8]  = PointPosition(29.183, 42.165, 16.864);
		atom_xyzs[9]  = PointPosition(30.136, 43.309, 16.704);
		atom_xyzs[10] = PointPosition(29.883, 44.187, 15.783);
		atom_xyzs[11] = PointPosition(28.978, 43.96, 14.678);
		atom_xyzs[12] = PointPosition(29.604, 43.507, 13.393);
		atom_xyzs[13] = PointPosition(30.563, 42.623, 13.495);
		atom_xyzs[14] = PointPosition(31.191, 42.012, 12.331);
		atom_xyzs[15] = PointPosition(30.459, 40.666, 12.13);

		jacobians[1] = 0.18665;
		jacobians[2] = 0.209767;

		conservative_solution = 2;
	}
};

class SixResidueTest : public TestHelper { // {{{1
public:
	SixResidueTest() : TestHelper("six", 6, 2, 3, 5, 4) {
		atom_xyzs[1]  = PointPosition(27.751, 35.867, 13.74);
		atom_xyzs[2]  = PointPosition(27.691, 37.315, 14.143);
		atom_xyzs[3]  = PointPosition(28.469, 37.475, 15.42);
		atom_xyzs[4]  = PointPosition(29.426, 38.43, 15.446);
		atom_xyzs[5]  = PointPosition(30.225, 38.643, 16.662);
		atom_xyzs[6]  = PointPosition(29.463, 38.096, 17.871);
		atom_xyzs[7]  = PointPosition(28.753, 38.941, 18.587);
		atom_xyzs[8]  = PointPosition(27.765, 39.841, 17.951);
		atom_xyzs[9]  = PointPosition(28.406, 40.842, 17.011);
		atom_xyzs[10] = PointPosition(29.023, 41.907, 17.556);
		atom_xyzs[11] = PointPosition(28.237, 43.188, 17.553);
		atom_xyzs[12] = PointPosition(27.353, 43.302, 16.35);
		atom_xyzs[13] = PointPosition(27.894, 43.145, 15.181);
		atom_xyzs[14] = PointPosition(28.978, 43.96, 14.678);
		atom_xyzs[15] = PointPosition(29.604, 43.507, 13.393);
		atom_xyzs[16] = PointPosition(30.563, 42.623, 13.495);
		atom_xyzs[17] = PointPosition(31.191, 42.012, 12.331);
		atom_xyzs[18] = PointPosition(30.459, 40.666, 12.13);

		jacobians[1] = 0.0623174;
		jacobians[2] = 0.0315628;
		jacobians[3] = 0.0392961;
		jacobians[4] = 0.0318445;

		conservative_solution = 1;
	}		
};

class SevenResidueTest : public TestHelper { // {{{1
public:
	SevenResidueTest() : TestHelper("seven", 7, 2, 4, 6, 10) {
		atom_xyzs[1]  = PointPosition(27.751, 35.867, 13.74);
		atom_xyzs[2]  = PointPosition(27.691, 37.315, 14.143);
		atom_xyzs[3]  = PointPosition(28.469, 37.475, 15.42);
		atom_xyzs[4]  = PointPosition(29.426, 38.43, 15.446);
		atom_xyzs[5]  = PointPosition(30.225, 38.643, 16.662);
		atom_xyzs[6]  = PointPosition(31.378, 37.637, 16.694);
		atom_xyzs[7]  = PointPosition(31.928, 37.292, 15.549);
		atom_xyzs[8]  = PointPosition(32.195, 38.296, 14.495);
		atom_xyzs[9]  = PointPosition(30.926, 38.899, 13.928);
		atom_xyzs[10] = PointPosition(29.817, 38.136, 13.924);
		atom_xyzs[11] = PointPosition(28.668, 38.65, 14.745);
		atom_xyzs[12] = PointPosition(28.082, 39.901, 14.168);
		atom_xyzs[13] = PointPosition(28.586, 40.356, 13.063);
		atom_xyzs[14] = PointPosition(28.214, 41.62, 12.466);
		atom_xyzs[15] = PointPosition(28.984, 42.821, 12.926);
		atom_xyzs[16] = PointPosition(30.142, 43.023, 12.352);
		atom_xyzs[17] = PointPosition(31.191, 42.012, 12.331);
		atom_xyzs[18] = PointPosition(30.459, 40.666, 12.13);
		atom_xyzs[19] = PointPosition(30.163, 40.338, 10.886);
		atom_xyzs[20] = PointPosition(29.542, 39.02, 10.653);
		atom_xyzs[21] = PointPosition(30.494, 38.261, 9.729);

		jacobians[1]  = 0.033083;
		jacobians[2]  = 0.0225511;
		jacobians[3]  = 0.0188586;
		jacobians[4]  = 0.0587687;
		jacobians[5]  = 0.120645;
		jacobians[6]  = 0.068554;
		jacobians[7]  = 0.0399261;
		jacobians[8]  = 0.0605512;
		jacobians[9]  = 0.0728305;
		jacobians[10] = 0.154974;

		conservative_solution = 7;
	}
};

class FoldTreeTest : public TestHelper { // {{{1
public:
	FoldTreeTest() : TestHelper("fold-tree", 8, 3, 4, 5, 2) {
		atom_xyzs[1]  = PointPosition(27.751, 35.867, 13.740);
		atom_xyzs[2]  = PointPosition(27.691, 37.315, 14.143);
		atom_xyzs[3]  = PointPosition(28.469, 37.475, 15.420);
		atom_xyzs[4]  = PointPosition(29.426, 38.430, 15.446);
		atom_xyzs[5]  = PointPosition(30.225, 38.643, 16.662);
		atom_xyzs[6]  = PointPosition(29.664, 39.839, 17.434);
		atom_xyzs[7]  = PointPosition(30.132, 40.069, 18.642);
		atom_xyzs[8]  = PointPosition(29.607, 41.180, 19.467);
		atom_xyzs[9]  = PointPosition(30.075, 42.538, 18.984);
		atom_xyzs[10] = PointPosition(30.991, 42.571, 17.998);
		atom_xyzs[11] = PointPosition(31.422, 43.940, 17.553);
		atom_xyzs[12] = PointPosition(30.755, 44.351, 16.277);
		atom_xyzs[13] = PointPosition(29.721, 43.673, 15.885);
		atom_xyzs[14] = PointPosition(28.978, 43.960, 14.678);
		atom_xyzs[15] = PointPosition(29.604, 43.507, 13.393);
		atom_xyzs[16] = PointPosition(30.563, 42.623, 13.495);
		atom_xyzs[17] = PointPosition(31.191, 42.012, 12.331);
		atom_xyzs[18] = PointPosition(30.459, 40.666, 12.130);
		atom_xyzs[19] = PointPosition(30.163, 40.338, 10.886);
		atom_xyzs[20] = PointPosition(29.542, 39.020, 10.653);
		atom_xyzs[21] = PointPosition(30.494, 38.261,  9.729);
		atom_xyzs[22] = PointPosition(30.795, 37.015, 10.095);
		atom_xyzs[23] = PointPosition(31.720, 36.289,  9.176);
		atom_xyzs[24] = PointPosition(30.955, 35.211,  8.459);

		jacobians[1]  = 0.2407;
		jacobians[2]  = 0.3203;

		conservative_solution = 2;
	}
};

class NumericalStabilityTest : public TestHelper { // {{{1
public:
	NumericalStabilityTest() : TestHelper("stability", 14, 2, 7, 12, 6) {
		atom_xyzs[1]  = PointPosition( 1.628,  0.206,  0.303);
		atom_xyzs[2]  = PointPosition( 2.217,  0.093, -1.029);
		atom_xyzs[3]  = PointPosition( 3.326, -0.933, -1.044);
		atom_xyzs[4]  = PointPosition( 3.758, -1.184, -2.331);
		atom_xyzs[5]  = PointPosition( 4.884, -2.108, -2.427);
		atom_xyzs[6]  = PointPosition( 5.426, -2.155, -3.835);
		atom_xyzs[7]  = PointPosition( 4.766, -1.377, -4.720);
		atom_xyzs[8]  = PointPosition( 5.201, -1.309, -6.112);
		atom_xyzs[9]  = PointPosition( 6.371, -0.367, -6.268);
		atom_xyzs[10] = PointPosition( 6.341,  0.322, -7.403);
		atom_xyzs[11] = PointPosition( 7.387,  1.296, -7.702);
		atom_xyzs[12] = PointPosition( 8.573,  1.115, -6.784);
		atom_xyzs[13] = PointPosition( 8.307,  0.607, -5.578);
		atom_xyzs[14] = PointPosition( 9.339,  0.345, -4.579);
		atom_xyzs[15] = PointPosition(10.378,  1.441, -4.576);
		atom_xyzs[16] = PointPosition(10.224,  2.331, -5.545);
		atom_xyzs[17] = PointPosition(11.168,  3.436, -5.691);
		atom_xyzs[18] = PointPosition(12.593,  2.941, -5.636);
		atom_xyzs[19] = PointPosition(13.065,  2.785, -4.393);
		atom_xyzs[20] = PointPosition(14.425,  2.307, -4.164);
		atom_xyzs[21] = PointPosition(14.711,  2.178, -2.688);
		atom_xyzs[22] = PointPosition(15.940,  1.758, -2.361);
		atom_xyzs[23] = PointPosition(16.382,  1.626, -0.976);
		atom_xyzs[24] = PointPosition(16.292,  0.191, -0.515);
		atom_xyzs[25] = PointPosition(15.072, -0.173, -0.082);
		atom_xyzs[26] = PointPosition(14.815, -1.549,  0.335);
		atom_xyzs[27] = PointPosition(13.619, -1.621,  1.253);
		atom_xyzs[28] = PointPosition(12.479, -1.862,  0.644);
		atom_xyzs[29] = PointPosition(11.245, -2.034,  1.405);
		atom_xyzs[30] = PointPosition(10.361, -0.817,  1.283);
		atom_xyzs[31] = PointPosition( 9.050, -1.043,  1.375);
		atom_xyzs[32] = PointPosition( 8.072,  0.038,  1.292);
		atom_xyzs[33] = PointPosition( 7.133,  0.010,  2.474);
		atom_xyzs[34] = PointPosition( 6.442, -1.089,  2.600);
		atom_xyzs[35] = PointPosition( 5.539, -1.306,  3.726);
		atom_xyzs[36] = PointPosition( 4.300, -0.452,  3.599);
		atom_xyzs[37] = PointPosition( 3.297, -0.624,  4.357);
		atom_xyzs[38] = PointPosition( 1.844, -0.766,  4.334);
		atom_xyzs[39] = PointPosition( 1.222,  0.201,  3.354);
		atom_xyzs[40] = PointPosition( 0.939,  1.145,  2.626);
		atom_xyzs[41] = PointPosition(-0.035,  1.500,  1.598);
		atom_xyzs[42] = PointPosition( 0.244,  0.760,  0.311);

		jacobians[1] = 0.00695385;
		jacobians[2] = 0.00704393;
		jacobians[3] = 0.00466329;
		jacobians[4] = 0.00427067;
		jacobians[5] = 0.0101787;
		jacobians[6] = 0.00797465;

		conservative_solution = 2;
	}
};
// }}}1

#endif
