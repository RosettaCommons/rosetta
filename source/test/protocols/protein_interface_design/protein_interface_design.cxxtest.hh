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
				// 	1-2 9
				// 	2-3 6
				// 	1-4 1
				// 	2-4 0
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

			void test_AtomicContactCount_filter()
			{
				core::pose::Pose testpose = fullatom_pose_from_string(AtomicContactCount_pose_string());

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
				protocols::toolbox::task_operations::PreventResiduesFromRepackingOperationOP prevent_3_repack( new protocols::toolbox::task_operations::PreventResiduesFromRepackingOperation(nopack_residues) );
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
			}
  };
}
