// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/protein_interface_design.cxxtest.hh
/// @brief Test suite for protocols/protein_interface_design/*
/// @author Alex Ford (fordas@uw.edu)
//
// Test headers
#include <iostream>
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <util/pose_funcs.hh>

// Project headers
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>


#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/task_operations/PreventResiduesFromRepackingOperation.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>

#include <protocols/protein_interface_design/filters/AtomicContactCountFilter.hh>

#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

static basic::Tracer TR("protocols.protein_interface_design.protein_interface_design.cxxtest");

namespace
{
/*
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::ResidueOP;
using core::conformation::ResidueFactory;
using core::chemical::ChemicalManager;
using core::chemical::ResidueType;
*/

class ProteinInterfaceDesignTests : public CxxTest::TestSuite {
public:
	void setUp()
	{
		core_init();
		// Initialize dock design parser to setup pose metric calculators.
		protocols::rosetta_scripts::RosettaScriptsParser parser;
	}

	std::string AtomicContactCount_pose_string()
	{
		// Test data for AtomicContactCounter filter
		// C-C <4.5 A contacts:
		//  1-2 9
		//  2-3 6
		//  1-4 1
		//  2-4 0
		return
			"ATOM      1  N   LEU A   1      60.238  41.546  48.182  1.00163.13           N  \n"
			"ATOM      2  CA  LEU A   1      58.985  41.560  47.367  1.00158.47           C  \n"
			"ATOM      3  C   LEU A   1      58.647  42.914  46.715  1.00155.76           C  \n"
			"ATOM      4  O   LEU A   1      57.823  43.683  47.216  1.00155.63           O  \n"
			"ATOM      5  CB  LEU A   1      57.778  41.005  48.133  1.00158.60           C  \n"
			"ATOM      6  CG  LEU A   1      56.462  40.871  47.368  1.00157.38           C  \n"
			"ATOM      7  CD1 LEU A   1      56.571  39.920  46.190  1.00156.08           C  \n"
			"ATOM      8  CD2 LEU A   1      55.407  40.398  48.319  1.00156.59           C  \n"
			"ATOM      9  H   LEU A   1      60.990  40.925  47.918  1.00  0.00           H  \n"
			"ATOM     10  HA  LEU A   1      59.122  40.948  46.476  1.00  0.00           H  \n"
			"ATOM     11 1HB  LEU A   1      58.159  40.013  48.369  1.00  0.00           H  \n"
			"ATOM     12 2HB  LEU A   1      57.599  41.555  49.057  1.00  0.00           H  \n"
			"ATOM     13  HG  LEU A   1      56.187  41.867  47.020  1.00  0.00           H  \n"
			"ATOM     14 1HD1 LEU A   1      55.609  39.860  45.680  1.00  0.00           H  \n"
			"ATOM     15 2HD1 LEU A   1      57.327  40.286  45.495  1.00  0.00           H  \n"
			"ATOM     16 3HD1 LEU A   1      56.854  38.930  46.547  1.00  0.00           H  \n"
			"ATOM     17 1HD2 LEU A   1      54.459  40.297  47.789  1.00  0.00           H  \n"
			"ATOM     18 2HD2 LEU A   1      55.696  39.432  48.733  1.00  0.00           H  \n"
			"ATOM     19 3HD2 LEU A   1      55.295  41.120  49.128  1.00  0.00           H  \n"
			"TER      20      LEU A   1\n"
			"ATOM     20  N   PHE B   2      50.173  39.998  52.534  1.00117.61           N  \n"
			"ATOM     21  CA  PHE B   2      51.538  40.388  52.863  1.00119.53           C  \n"
			"ATOM     22  C   PHE B   2      51.652  41.811  53.394  1.00120.68           C  \n"
			"ATOM     23  O   PHE B   2      52.124  42.010  54.512  1.00120.20           O  \n"
			"ATOM     24  CB  PHE B   2      52.410  40.199  51.630  1.00119.38           C  \n"
			"ATOM     25  CG  PHE B   2      53.803  40.769  51.751  1.00120.47           C  \n"
			"ATOM     26  CD1 PHE B   2      54.049  42.122  51.500  1.00120.19           C  \n"
			"ATOM     27  CD2 PHE B   2      54.884  39.938  52.055  1.00120.92           C  \n"
			"ATOM     28  CE1 PHE B   2      55.340  42.640  51.577  1.00120.28           C  \n"
			"ATOM     29  CE2 PHE B   2      56.182  40.441  52.136  1.00120.98           C  \n"
			"ATOM     30  CZ  PHE B   2      56.410  41.792  51.897  1.00120.50           C  \n"
			"ATOM     31  H   PHE B   2      49.925  39.842  51.567  1.00  0.00           H  \n"
			"ATOM     32  HA  PHE B   2      51.915  39.765  53.675  1.00  0.00           H  \n"
			"ATOM     33 1HB  PHE B   2      52.534  39.138  51.417  1.00  0.00           H  \n"
			"ATOM     34 2HB  PHE B   2      51.955  40.690  50.771  1.00  0.00           H  \n"
			"ATOM     35  HD1 PHE B   2      53.216  42.776  51.242  1.00  0.00           H  \n"
			"ATOM     36  HD2 PHE B   2      54.697  38.880  52.241  1.00  0.00           H  \n"
			"ATOM     37  HE1 PHE B   2      55.514  43.699  51.389  1.00  0.00           H  \n"
			"ATOM     38  HE2 PHE B   2      57.015  39.784  52.384  1.00  0.00           H  \n"
			"ATOM     39  HZ  PHE B   2      57.422  42.189  51.958  1.00  0.00           H  \n"
			"ATOM     56  N   SER B   3      55.588  44.222  55.854  1.00  0.00           N  \n"
			"ATOM     57  CA  SER B   3      55.371  42.778  55.825  1.00  0.00           C  \n"
			"ATOM     58  C   SER B   3      56.155  42.091  56.918  1.00  0.00           C  \n"
			"ATOM     59  O   SER B   3      56.867  42.719  57.703  1.00  0.00           O  \n"
			"ATOM     60  CB  SER B   3      55.712  42.198  54.429  1.00  0.00           C  \n"
			"ATOM     61  OG  SER B   3      57.113  42.231  54.139  1.00  0.00           O  \n"
			"ATOM     62  H   SER B   3      56.233  44.690  56.582  1.00  0.00           H  \n"
			"ATOM     63  HA  SER B   3      54.302  42.589  56.031  1.00  0.00           H  \n"
			"ATOM     64 2HB  SER B   3      55.157  42.725  53.630  1.00  0.00           H  \n"
			"ATOM     65 3HB  SER B   3      55.368  41.147  54.369  1.00  0.00           H  \n"
			"ATOM     66  HG  SER B   3      57.381  43.153  54.123  1.00  0.00           H  \n"
			"TER      66      SER B   3\n"
			"ATOM     40  N   VAL C   4      55.443  38.449  39.966  1.00  0.00           N  \n"
			"ATOM     41  CA  VAL C   4      54.313  37.732  40.551  1.00  0.00           C  \n"
			"ATOM     42  C   VAL C   4      51.940  38.344  40.623  1.00  0.00           C  \n"
			"ATOM     43  O   VAL C   4      52.308  37.679  39.647  1.00  0.00           O  \n"
			"ATOM     44  CB  VAL C   4      53.907  38.401  41.924  1.00  0.00           C  \n"
			"ATOM     45  CG1 VAL C   4      52.605  37.875  42.589  1.00  0.00           C  \n"
			"ATOM     46  CG2 VAL C   4      54.995  38.295  43.019  1.00  0.00           C  \n"
			"ATOM     47  H   VAL C   4      55.395  38.884  38.978  1.00  0.00           H  \n"
			"ATOM     48  HA  VAL C   4      54.612  36.684  40.730  1.00  0.00           H  \n"
			"ATOM     49  HB  VAL C   4      53.745  39.483  41.727  1.00  0.00           H  \n"
			"ATOM     50 1HG1 VAL C   4      52.651  36.792  42.811  1.00  0.00           H  \n"
			"ATOM     51 2HG1 VAL C   4      52.378  38.396  43.539  1.00  0.00           H  \n"
			"ATOM     52 3HG1 VAL C   4      51.712  38.039  41.955  1.00  0.00           H  \n"
			"ATOM     53 1HG2 VAL C   4      55.963  38.707  42.679  1.00  0.00           H  \n"
			"ATOM     54 2HG2 VAL C   4      54.723  38.868  43.927  1.00  0.00           H  \n"
			"ATOM     55 3HG2 VAL C   4      55.662  37.430  43.197  1.00  0.00           H  \n"
			"END\n";
	}

	std::string AtomicContactCount_pose2_string()
	{
		return
			"ATOM      1  N   LYS A   6       1.272  21.035   1.792  1.00  0.00           N  \n"
			"ATOM      2  CA  LYS A   6       0.648  20.378   0.652  1.00  0.00           C  \n"
			"ATOM      3  C   LYS A   6       1.165  18.960   0.479  1.00  0.00           C  \n"
			"ATOM      4  O   LYS A   6       0.386  18.047   0.201  1.00  0.00           O  \n"
			"ATOM      5  CB  LYS A   6       0.868  21.184  -0.627  1.00  0.00           C  \n"
			"ATOM      6  CG  LYS A   6       0.056  22.480  -0.700  1.00  0.00           C  \n"
			"ATOM      7  CD  LYS A   6       0.346  23.244  -1.983  1.00  0.00           C  \n"
			"ATOM      8  CE  LYS A   6      -0.438  24.543  -2.063  1.00  0.00           C  \n"
			"ATOM      9  NZ  LYS A   6      -0.161  25.274  -3.335  1.00  0.00           N1+\n"
			"ATOM     10  H   LYS A   6       1.825  21.878   1.651  1.00  0.00           H  \n"
			"ATOM     11  HA  LYS A   6      -0.426  20.319   0.839  1.00  0.00           H  \n"
			"ATOM     12 1HB  LYS A   6       1.926  21.450  -0.706  1.00  0.00           H  \n"
			"ATOM     13 2HB  LYS A   6       0.617  20.572  -1.493  1.00  0.00           H  \n"
			"ATOM     14 1HG  LYS A   6      -1.008  22.239  -0.662  1.00  0.00           H  \n"
			"ATOM     15 2HG  LYS A   6       0.293  23.108   0.155  1.00  0.00           H  \n"
			"ATOM     16 1HD  LYS A   6       1.411  23.482  -2.019  1.00  0.00           H  \n"
			"ATOM     17 2HD  LYS A   6       0.098  22.624  -2.845  1.00  0.00           H  \n"
			"ATOM     18 1HE  LYS A   6      -1.503  24.327  -2.001  1.00  0.00           H  \n"
			"ATOM     19 2HE  LYS A   6      -0.156  25.179  -1.222  1.00  0.00           H  \n"
			"ATOM     20 1HZ  LYS A   6      -0.691  26.133  -3.356  1.00  0.00           H  \n"
			"ATOM     21 2HZ  LYS A   6       0.828  25.487  -3.395  1.00  0.00           H  \n"
			"ATOM     22 3HZ  LYS A   6      -0.431  24.695  -4.121  1.00  0.00           H  \n"
			"ATOM     23  N   ALA A   7       2.471  18.759   0.656  1.00  0.00           N  \n"
			"ATOM     24  CA  ALA A   7       3.036  17.426   0.504  1.00  0.00           C  \n"
			"ATOM     25  C   ALA A   7       2.433  16.459   1.519  1.00  0.00           C  \n"
			"ATOM     26  O   ALA A   7       2.102  15.318   1.183  1.00  0.00           O  \n"
			"ATOM     27  CB  ALA A   7       4.540  17.487   0.667  1.00  0.00           C  \n"
			"ATOM     28  H   ALA A   7       3.084  19.546   0.865  1.00  0.00           H  \n"
			"ATOM     29  HA  ALA A   7       2.797  17.062  -0.494  1.00  0.00           H  \n"
			"ATOM     30 1HB  ALA A   7       4.960  16.504   0.536  1.00  0.00           H  \n"
			"ATOM     31 2HB  ALA A   7       4.955  18.165  -0.079  1.00  0.00           H  \n"
			"ATOM     32 3HB  ALA A   7       4.776  17.853   1.666  1.00  0.00           H  \n"
			"ATOM     33  N   MET A   8       2.237  16.926   2.755  1.00  0.00           N  \n"
			"ATOM     34  CA  MET A   8       1.643  16.070   3.774  1.00  0.00           C  \n"
			"ATOM     35  C   MET A   8       0.210  15.704   3.414  1.00  0.00           C  \n"
			"ATOM     36  O   MET A   8      -0.194  14.544   3.544  1.00  0.00           O  \n"
			"ATOM     37  CB  MET A   8       1.678  16.755   5.135  1.00  0.00           C  \n"
			"ATOM     38  CG  MET A   8       3.064  16.875   5.744  1.00  0.00           C  \n"
			"ATOM     39  SD  MET A   8       3.029  17.426   7.458  1.00  0.00           S  \n"
			"ATOM     40  CE  MET A   8       2.611  19.158   7.265  1.00  0.00           C  \n"
			"ATOM     41  H   MET A   8       2.544  17.867   2.994  1.00  0.00           H  \n"
			"ATOM     42  HA  MET A   8       2.216  15.148   3.828  1.00  0.00           H  \n"
			"ATOM     43 1HB  MET A   8       1.273  17.759   5.031  1.00  0.00           H  \n"
			"ATOM     44 2HB  MET A   8       1.040  16.217   5.834  1.00  0.00           H  \n"
			"ATOM     45 1HG  MET A   8       3.580  15.924   5.695  1.00  0.00           H  \n"
			"ATOM     46 2HG  MET A   8       3.640  17.589   5.171  1.00  0.00           H  \n"
			"ATOM     47 1HE  MET A   8       2.567  19.632   8.245  1.00  0.00           H  \n"
			"ATOM     48 2HE  MET A   8       3.368  19.654   6.655  1.00  0.00           H  \n"
			"ATOM     49 3HE  MET A   8       1.641  19.244   6.784  1.00  0.00           H  \n"
			"TER      50      MET A   8\n"
			"ATOM     51  N   VAL B  68       8.058  16.020   7.969  1.00  0.00           N  \n"
			"ATOM     52  CA  VAL B  68       8.142  16.549   6.612  1.00  0.00           C  \n"
			"ATOM     53  C   VAL B  68       9.288  17.549   6.547  1.00  0.00           C  \n"
			"ATOM     54  O   VAL B  68       9.487  18.310   7.495  1.00  0.00           O  \n"
			"ATOM     55  CB  VAL B  68       6.807  17.224   6.205  1.00  0.00           C  \n"
			"ATOM     56  CG1 VAL B  68       6.495  18.391   7.136  1.00  0.00           C  \n"
			"ATOM     57  CG2 VAL B  68       6.878  17.722   4.753  1.00  0.00           C  \n"
			"ATOM     58  H   VAL B  68       8.185  16.657   8.743  1.00  0.00           H  \n"
			"ATOM     59  HA  VAL B  68       8.349  15.730   5.923  1.00  0.00           H  \n"
			"ATOM     60  HB  VAL B  68       6.006  16.494   6.300  1.00  0.00           H  \n"
			"ATOM     61 1HG1 VAL B  68       5.552  18.837   6.847  1.00  0.00           H  \n"
			"ATOM     62 2HG1 VAL B  68       6.426  18.034   8.160  1.00  0.00           H  \n"
			"ATOM     63 3HG1 VAL B  68       7.274  19.148   7.067  1.00  0.00           H  \n"
			"ATOM     64 1HG2 VAL B  68       5.927  18.175   4.481  1.00  0.00           H  \n"
			"ATOM     65 2HG2 VAL B  68       7.666  18.463   4.648  1.00  0.00           H  \n"
			"ATOM     66 3HG2 VAL B  68       7.080  16.888   4.089  1.00  0.00           H  \n"
			"TER      67      VAL B  68\n"
			"ATOM     68  N   TRP C 109       1.904  20.372  -6.341  1.00  0.00           N  \n"
			"ATOM     69  CA  TRP C 109       2.690  19.240  -5.863  1.00  0.00           C  \n"
			"ATOM     70  C   TRP C 109       2.041  17.903  -6.189  1.00  0.00           C  \n"
			"ATOM     71  O   TRP C 109       0.817  17.754  -6.118  1.00  0.00           O  \n"
			"ATOM     72  CB  TRP C 109       2.895  19.326  -4.348  1.00  0.00           C  \n"
			"ATOM     73  CG  TRP C 109       3.822  20.400  -3.930  1.00  0.00           C  \n"
			"ATOM     74  CD1 TRP C 109       3.478  21.637  -3.508  1.00  0.00           C  \n"
			"ATOM     75  CD2 TRP C 109       5.265  20.350  -3.890  1.00  0.00           C  \n"
			"ATOM     76  CE2 TRP C 109       5.701  21.579  -3.442  1.00  0.00           C  \n"
			"ATOM     77  CE3 TRP C 109       6.198  19.385  -4.206  1.00  0.00           C  \n"
			"ATOM     78  NE1 TRP C 109       4.594  22.359  -3.213  1.00  0.00           N  \n"
			"ATOM     79  CZ2 TRP C 109       7.050  21.867  -3.299  1.00  0.00           C  \n"
			"ATOM     80  CZ3 TRP C 109       7.543  19.685  -4.064  1.00  0.00           C  \n"
			"ATOM     81  CH2 TRP C 109       7.950  20.887  -3.627  1.00  0.00           C  \n"
			"ATOM     82  H   TRP C 109       1.278  20.839  -5.700  1.00  0.00           H  \n"
			"ATOM     83  HA  TRP C 109       3.664  19.271  -6.344  1.00  0.00           H  \n"
			"ATOM     84 1HB  TRP C 109       1.939  19.504  -3.864  1.00  0.00           H  \n"
			"ATOM     85 2HB  TRP C 109       3.280  18.376  -3.978  1.00  0.00           H  \n"
			"ATOM     86  HD1 TRP C 109       2.464  22.002  -3.425  1.00  0.00           H  \n"
			"ATOM     87  HE1 TRP C 109       4.600  23.335  -2.882  1.00  0.00           H  \n"
			"ATOM     88  HE3 TRP C 109       5.882  18.424  -4.571  1.00  0.00           H  \n"
			"ATOM     89  HZ2 TRP C 109       7.411  22.829  -2.957  1.00  0.00           H  \n"
			"ATOM     90  HZ3 TRP C 109       8.269  18.942  -4.318  1.00  0.00           H  \n"
			"ATOM     91  HH2 TRP C 109       9.019  21.086  -3.535  1.00  0.00           H  \n"
			"END\n";
	}

	void test_AtomicContactCount_filter()
	{
		core::pose::Pose testpose = fullatom_pose_from_string(AtomicContactCount_pose_string());
		core::pose::Pose testpose2 = fullatom_pose_from_string(AtomicContactCount_pose2_string());

		using namespace protocols::protein_interface_design::filters;

		// All contacts in pose
		AtomicContactCountFilter default_filter;
		core::Real default_count = default_filter.compute(testpose);
		TS_ASSERT_EQUALS(default_count, 16);

		// All cross-chain contacts in pose
		AtomicContactCountFilter cross_chain_filter;
		cross_chain_filter.initialize_cross_chain();
		core::Real cross_chain_count = cross_chain_filter.compute(testpose);
		TS_ASSERT_EQUALS(cross_chain_count, 10);

		// Test jump specification
		// By default, fold tree connects jump 1 between res 1&2 and jump 2 between res 1&4
		AtomicContactCountFilter jump1_filter;
		jump1_filter.initialize_cross_jump(1);
		core::Real jump1_count = jump1_filter.compute(testpose);
		TS_ASSERT_EQUALS(jump1_count, 9);

		AtomicContactCountFilter jump2_filter;
		jump2_filter.initialize_cross_jump(2);
		core::Real jump2_count = jump2_filter.compute(testpose);
		TS_ASSERT_EQUALS(jump2_count, 1);

		// Test with TaskOperation restriction
		core::pack::task::TaskFactoryOP select_12( new core::pack::task::TaskFactory() );
		utility::vector1<core::Size> nopack_residues;
		nopack_residues.push_back(3);
		nopack_residues.push_back(4);
		protocols::task_operations::PreventResiduesFromRepackingOperationOP prevent_3_repack( new protocols::task_operations::PreventResiduesFromRepackingOperation(nopack_residues) );
		select_12->push_back(prevent_3_repack);

		AtomicContactCountFilter res12_filter;
		res12_filter.initialize_all_atoms(select_12);
		core::Real res12_count = res12_filter.compute(testpose);
		TS_ASSERT_EQUALS(res12_count, 9);

		// Test with SASA normalization, partition by jump
		protocols::simple_filters::InterfaceSasaFilter sasafilter;
		sasafilter.jump(1);
		core::Real jump1_sasa = sasafilter.compute(testpose);

		protocols::simple_filters::InterfaceSasaFilter sasafilter2;
		sasafilter2.jump(2);
		core::Real jump2_sasa = sasafilter2.compute(testpose);

		AtomicContactCountFilter sasanorm_filter;
		sasanorm_filter.initialize_cross_jump(1, "", NULL, true);
		core::Real sasanorm_count = sasanorm_filter.compute(testpose);
		TS_ASSERT_DELTA(9 / jump1_sasa, sasanorm_count, 1e-6);

		// Test with SASA normalization, partition all chains
		AtomicContactCountFilter sasanorm_filter_chains;
		sasanorm_filter_chains.initialize_cross_chain(NULL, true);
		core::Real sasanorm_chains_count = sasanorm_filter_chains.compute(testpose);
		TS_ASSERT_DELTA(10 / (jump1_sasa + jump2_sasa), sasanorm_chains_count, 1e-6);

		// Test with SASA normalization, partition all chains limited to task operation
		AtomicContactCountFilter sasanorm_filter_chains_selection;
		sasanorm_filter_chains_selection.initialize_cross_chain(select_12, true);
		core::Real sasanorm_chains_selection_count = sasanorm_filter_chains_selection.compute(testpose);
		TS_ASSERT_DELTA(9 / (jump1_sasa + jump2_sasa), sasanorm_chains_selection_count, 1e-6);

		// Test with SASA normalization, partition all chains limited to task operation and detect chain separation
		AtomicContactCountFilter sasanorm_filter_chains_selection_detect;
		sasanorm_filter_chains_selection_detect.initialize_cross_chain(select_12, true, true);
		core::Real sasanorm_chains_selection_detect_count = sasanorm_filter_chains_selection_detect.compute(testpose);
		TS_ASSERT_DELTA(9 / (jump1_sasa), sasanorm_chains_selection_detect_count, 1e-6);

		// Test without contact
		// Separate jump 2
		protocols::rigid::RigidBodyTransMover jump2_translate( testpose, 2, false );
		jump2_translate.step_size( 1000.0 );
		jump2_translate.trans_axis(core::Vector(1, 0, 0));
		jump2_translate.apply( testpose );

		core::Real default_count_sep = default_filter.compute(testpose);
		TS_ASSERT_EQUALS(default_count_sep, 15);

		core::Real jump2_count_sep = jump2_filter.compute(testpose);
		TS_ASSERT_EQUALS(jump2_count_sep, 0);

		core::Real cross_chain_count_sep = cross_chain_filter.compute(testpose);
		TS_ASSERT_EQUALS(cross_chain_count_sep, 9);

		// Separate jump 1
		protocols::rigid::RigidBodyTransMover jump1_translate( testpose, 1, false );
		jump1_translate.step_size( 1000.0 );
		jump1_translate.trans_axis(core::Vector(1, 0, 0));
		jump1_translate.apply( testpose );

		core::Real default_count_allsep = default_filter.compute(testpose);
		TS_ASSERT_EQUALS(default_count_allsep, 6);

		core::Real sasanorm_count_allsep = sasanorm_filter.compute(testpose);
		TS_ASSERT_EQUALS(0, sasanorm_count_allsep);

		core::Real jump2_count_allsep = jump2_filter.compute(testpose);
		TS_ASSERT_EQUALS(jump2_count_allsep, 0);

		core::Real cross_chain_count_allsep = cross_chain_filter.compute(testpose);
		TS_ASSERT_EQUALS(cross_chain_count_allsep, 0);

		// Test non_local, res_contact and count_SD_NE1 options:
		// One contact per residue, local, counting SD and NE1 atoms:
		AtomicContactCountFilter non_local_res_filter;
		non_local_res_filter.initialize_all_atoms( NULL, false, NULL, false, false, true, true);
		core::Real non_local_res_filter_count = non_local_res_filter.compute(testpose2);
		TS_ASSERT_EQUALS( non_local_res_filter_count, 2 );

		// Counting all contacts per residue, local, counting SD and NE1 atoms:
		AtomicContactCountFilter hphob_sc_contacts_local;
		hphob_sc_contacts_local.initialize_all_atoms( NULL, false, NULL, false, false, false, true);
		core::Real hphob_sc_contacts_local_count = hphob_sc_contacts_local.compute(testpose2);
		TS_ASSERT_EQUALS( hphob_sc_contacts_local_count, 8 );

		// Counting all contacts per residue, non local, not counting SD and NE1 atoms:
		AtomicContactCountFilter hphob_sc_contacts_no_N_S_nonlocal;
		hphob_sc_contacts_no_N_S_nonlocal.initialize_all_atoms( NULL, false, NULL, false, true, false, false);
		core::Real hphob_sc_contacts_no_N_S_nonlocal_count = hphob_sc_contacts_no_N_S_nonlocal.compute(testpose2);
		TS_ASSERT_EQUALS( hphob_sc_contacts_no_N_S_nonlocal_count, 6 );

	}
};
}
