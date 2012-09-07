// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>

#include <protocols/protein_interface_design/filters/AtomicContactCountFilter.hh>

#include <protocols/jd2/DockDesignParser.hh>

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
				protocols::jd2::DockDesignParser parser;
      }

			std::string AtomicContactCount_pose_string()
			{
				// Test data for AtomicContactCounter filter
				// C-C <4.5 A contacts:
				// 	1-2 9
				// 	1-3 1
				// 	2-3 0
				return 
				"ATOM      1  N   LEU A   1      60.449  43.365  47.194  1.00163.13           N\n"
					"ATOM      2  CA  LEU A   1      59.196  43.379  46.379  1.00158.47           C\n"
					"ATOM      3  C   LEU A   1      58.858  44.733  45.727  1.00155.76           C\n"
					"ATOM      4  O   LEU A   1      58.034  45.502  46.228  1.00155.63           O\n"
					"ATOM      5  CB  LEU A   1      57.989  42.824  47.145  1.00158.60           C\n"
					"ATOM      6  CG  LEU A   1      56.673  42.690  46.380  1.00157.38           C\n"
					"ATOM      7  CD1 LEU A   1      56.782  41.739  45.202  1.00156.08           C\n"
					"ATOM      8  CD2 LEU A   1      55.618  42.217  47.331  1.00156.59           C\n"
					"ATOM      9  H   LEU A   1      61.201  42.744  46.930  1.00  0.00           H\n"
					"ATOM     10  HA  LEU A   1      59.333  42.767  45.488  1.00  0.00           H\n"
					"ATOM     11 1HB  LEU A   1      58.370  41.832  47.381  1.00  0.00           H\n"
					"ATOM     12 2HB  LEU A   1      57.810  43.374  48.069  1.00  0.00           H\n"
					"ATOM     13  HG  LEU A   1      56.398  43.686  46.032  1.00  0.00           H\n"
					"ATOM     14 1HD1 LEU A   1      55.820  41.679  44.692  1.00  0.00           H\n"
					"ATOM     15 2HD1 LEU A   1      57.538  42.105  44.507  1.00  0.00           H\n"
					"ATOM     16 3HD1 LEU A   1      57.065  40.749  45.559  1.00  0.00           H\n"
					"ATOM     17 1HD2 LEU A   1      54.670  42.116  46.801  1.00  0.00           H\n"
					"ATOM     18 2HD2 LEU A   1      55.907  41.251  47.745  1.00  0.00           H\n"
					"ATOM     19 3HD2 LEU A   1      55.506  42.939  48.140  1.00  0.00           H\n"
					"TER      20      LEU A   1\n"
					"ATOM     21  N   PHE B   2      50.384  41.817  51.546  1.00117.61           N\n"
					"ATOM     22  CA  PHE B   2      51.749  42.207  51.875  1.00119.53           C\n"
					"ATOM     23  C   PHE B   2      51.863  43.630  52.406  1.00120.68           C\n"
					"ATOM     24  O   PHE B   2      52.335  43.829  53.524  1.00120.20           O\n"
					"ATOM     25  CB  PHE B   2      52.621  42.018  50.642  1.00119.38           C\n"
					"ATOM     26  CG  PHE B   2      54.014  42.588  50.763  1.00120.47           C\n"
					"ATOM     27  CD1 PHE B   2      54.260  43.941  50.512  1.00120.19           C\n"
					"ATOM     28  CD2 PHE B   2      55.095  41.757  51.067  1.00120.92           C\n"
					"ATOM     29  CE1 PHE B   2      55.551  44.459  50.589  1.00120.28           C\n"
					"ATOM     30  CE2 PHE B   2      56.393  42.260  51.148  1.00120.98           C\n"
					"ATOM     31  CZ  PHE B   2      56.621  43.611  50.909  1.00120.50           C\n"
					"ATOM     32  H   PHE B   2      50.136  41.661  50.579  1.00  0.00           H\n"
					"ATOM     33  HA  PHE B   2      52.126  41.584  52.687  1.00  0.00           H\n"
					"ATOM     34 1HB  PHE B   2      52.745  40.957  50.429  1.00  0.00           H\n"
					"ATOM     35 2HB  PHE B   2      52.166  42.509  49.783  1.00  0.00           H\n"
					"ATOM     36  HD1 PHE B   2      53.427  44.595  50.254  1.00  0.00           H\n"
					"ATOM     37  HD2 PHE B   2      54.908  40.699  51.253  1.00  0.00           H\n"
					"ATOM     38  HE1 PHE B   2      55.725  45.518  50.401  1.00  0.00           H\n"
					"ATOM     39  HE2 PHE B   2      57.226  41.603  51.396  1.00  0.00           H\n"
					"ATOM     40  HZ  PHE B   2      57.633  44.008  50.970  1.00  0.00           H\n"
					"TER      41      PHE B   2\n"
					"ATOM     42  N   VAL C   3      55.654  40.268  38.978  1.00  0.00           N\n"
					"ATOM     43  CA  VAL C   3      54.524  39.551  39.563  1.00  0.00           C\n"
					"ATOM     44  C   VAL C   3      52.151  40.163  39.635  1.00  0.00           C\n"
					"ATOM     45  O   VAL C   3      52.519  39.498  38.659  1.00  0.00           O\n"
					"ATOM     46  CB  VAL C   3      54.118  40.220  40.936  1.00  0.00           C\n"
					"ATOM     47  CG1 VAL C   3      52.816  39.694  41.601  1.00  0.00           C\n"
					"ATOM     48  CG2 VAL C   3      55.206  40.114  42.031  1.00  0.00           C\n"
					"ATOM     49  H   VAL C   3      55.606  40.703  37.990  1.00  0.00           H\n"
					"ATOM     50  HA  VAL C   3      54.823  38.503  39.742  1.00  0.00           H\n"
					"ATOM     51  HB  VAL C   3      53.956  41.302  40.739  1.00  0.00           H\n"
					"ATOM     52 1HG1 VAL C   3      52.862  38.611  41.823  1.00  0.00           H\n"
					"ATOM     53 2HG1 VAL C   3      52.589  40.215  42.551  1.00  0.00           H\n"
					"ATOM     54 3HG1 VAL C   3      51.923  39.858  40.967  1.00  0.00           H\n"
					"ATOM     55 1HG2 VAL C   3      56.174  40.526  41.691  1.00  0.00           H\n"
					"ATOM     56 2HG2 VAL C   3      54.934  40.687  42.939  1.00  0.00           H\n"
					"ATOM     57 3HG2 VAL C   3      55.873  39.249  42.209  1.00  0.00           H\n"
					"END\n";
			}

			void test_AtomicContactCount_filter()
			{
				core::pose::Pose testpose = fullatom_pose_from_string(AtomicContactCount_pose_string());

				// Rebuild pose fold tree connectivity.
				core::kinematics::FoldTree ft = testpose.fold_tree();
				ft.clear();
				ft.add_edge(1, 2, 1);
				ft.add_edge(2, 3, 2);
				ft.reorder( 1 );
				testpose.fold_tree( ft );

				using namespace protocols::protein_interface_design::filters;

				// Test default case, all contacts across jump 1
				AtomicContactCountFilter default_filter;
				core::Real default_count = default_filter.compute(testpose);
				TS_ASSERT_EQUALS(default_count, 10);

				// Test alternate jump
				AtomicContactCountFilter jump2_filter(NULL, 2, 4.5, false);
				core::Real jump2_count = jump2_filter.compute(testpose);
				TS_ASSERT_EQUALS(jump2_count, 1);

				// Test with TaskOperation restriction
				core::pack::task::TaskFactoryOP select_12 = new core::pack::task::TaskFactory();
				utility::vector1<core::Size> nopack_residues;
				nopack_residues.push_back(3);
				protocols::toolbox::task_operations::PreventResiduesFromRepackingOperationOP prevent_3_repack = 
					new protocols::toolbox::task_operations::PreventResiduesFromRepackingOperation(nopack_residues);
				select_12->push_back(prevent_3_repack);

				AtomicContactCountFilter res12_filter(select_12, 1, 4.5, false);
				core::Real res12_count = res12_filter.compute(testpose);
				TS_ASSERT_EQUALS(res12_count, 9);

				// Test with SASA normalization
				protocols::simple_filters::InterfaceSasaFilter sasafilter;
				sasafilter.jump(1);
				core::Real jump1_sasa = sasafilter.compute(testpose);

				AtomicContactCountFilter sasanorm_filter(NULL, 1, 4.5, true);
				core::Real sasanorm_count = sasanorm_filter.compute(testpose);
				TS_ASSERT_EQUALS(10 / jump1_sasa, sasanorm_count);

				// Test without contact
				// Separate jump 2
				protocols::rigid::RigidBodyTransMover jump2_translate( testpose, 2);
				jump2_translate.step_size( 1000.0 );
				jump2_translate.trans_axis(core::Vector(1, 0, 0));
				jump2_translate.apply( testpose );
				
				core::Real default_count_sep = default_filter.compute(testpose);
				TS_ASSERT_EQUALS(default_count_sep, 9);

				core::Real jump2_count_sep = jump2_filter.compute(testpose);
				TS_ASSERT_EQUALS(jump2_count_sep, 0);

				// Separate jump 1
				protocols::rigid::RigidBodyTransMover jump1_translate( testpose, 1);
				jump1_translate.step_size( 1000.0 );
				jump1_translate.trans_axis(core::Vector(1, 0, 0));
				jump1_translate.apply( testpose );

				core::Real default_count_allsep = default_filter.compute(testpose);
				TS_ASSERT_EQUALS(default_count_allsep, 0);

				core::Real sasanorm_count_allsep = sasanorm_filter.compute(testpose);
				TS_ASSERT_EQUALS(0, sasanorm_count_allsep);

				core::Real jump2_count_allsep = jump2_filter.compute(testpose);
				TS_ASSERT_EQUALS(jump2_count_allsep, 0);
			}
  };
}
