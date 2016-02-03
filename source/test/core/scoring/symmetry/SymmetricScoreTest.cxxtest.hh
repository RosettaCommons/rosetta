// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/symmetry/ScoreTest.cxxtest.hh
/// @brief  unified scoring test.
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/scoring/Energies.hh>
#include <core/types.hh>

// Unit headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/symmetry/util.hh>

// Package headers

#include <test/UTracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <utility/py/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <execinfo.h>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.symmetry.SymmetricScoreTest.cxxtest");

// using declarations
using namespace core;
using namespace scoring;
using namespace scoring::symmetry;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;

///////////////////////////////////////////////////////////////////////////
/// @name ScoreTest
/// @brief: unified tests for difference score functions/methods
///////////////////////////////////////////////////////////////////////////
class SymmetricScoreTest : public CxxTest::TestSuite {

public:
	chemical::ResidueTypeSetCAP residue_set;

	void setUp() {
		core_init_with_additional_options( "-no_optH -symmetry:symmetry_definition core/scoring/symmetry/sym_def.dat" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	void tearDown() {}

	void test_SymmetricScoreTest() {
		//one_score_type_test(, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/.u");

		one_score_type_test(fa_atr, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_atr.data");

		one_score_type_test(fa_rep, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_rep.data");
		one_score_type_test(fa_sol, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_sol.data");

		one_score_type_test(fa_intra_atr, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_intra_atr.data");
		one_score_type_test(fa_intra_rep, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_intra_rep.data");
		one_score_type_test(fa_intra_sol, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_intra_sol.data");

		one_score_type_test(mm_twist, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/mm_twist.data");
		//one_score_type_test(mm_bend, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/mm_bend.u");
		//one_score_type_test(mm_stretch, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/mm_stretch.u");
		one_score_type_test(fa_elec, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_elec.data");
		one_score_type_test(atom_pair_constraint, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/atom_pair_constraint.data");
		one_score_type_test(coordinate_constraint, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/coordinate_constraint.data");

		// 0 one_score_type_test(angle_constraint, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/angle_constraint.u");
		// 0 one_score_type_test(dihedral_constraint, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/dihedral_constraint.u");
		// 0 one_score_type_test(dna_bp, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/dna_bp.u");
		// 0 one_score_type_test(dna_bs, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/dna_bs.u");

		// nw one_score_type_test(fa_pair, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_pair.u");
		//one_score_type_test(fa_plane, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/fa_plane.u");
		one_score_type_test(hbond_sr_bb, "core/scoring/symmetry/test_in.pdb", "core/scoring/symmetry/hbond_sr_bb.data");
		/*
		hbond_lr_bb,
		hbond_bb_sc,
		hbond_sc,
		gb_elec,
		dslf_ss_dst,
		dslf_cs_ang,
		dslf_ss_dih,
		dslf_ca_dih,
		dslf_cbs_ds,

		rama,
		omega,
		fa_dun,
		p_aa_pp,
		ref,
		envsmooth,

		// rigid body move specific scores begin - Monica Berrondo
		rb_scorefxn,
		rb_env,
		rb_pair,
		rb_cont,
		rb_cont_cap,
		rb_vdw,
		rb_site_cst,
		rb_fab,
		rb_fab_cap,
		rb_wsl_elec,
		// rigid body move specific scores end - Monica Berrondo

		// centroid scores
		env,
		pair,
		cbeta,
		vdw,
		rg,
		cenpack,
		hs_pair,
		ss_pair,
		rsigma,


		chainbreak,
		*/
		// test score difference
		test_sym_vs_asym_score();
	}

	void one_score_type_test(scoring::ScoreType st, std::string pdb_file_name, std::string data_file_name,
		double abs_p=0.0001, double rel_p=0.0001)
	{
		//std::cout << " Testing score: " << scoring::name_from_score_type(st) << "..." << std::endl;

		Pose pose;
		core::import_pose::pose_from_file(pose, pdb_file_name, core::import_pose::PDB_file);
		core::pose::symmetry::make_symmetric_pose( pose );

		SymmetricScoreFunction scorefxn;
		scorefxn.set_weight(st, 1.0 );

		Energy score = scorefxn( pose );
		/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( pose );

		std::vector<double> D;

		D.push_back( score );

		for ( Size r=1; r<=pose.total_residue(); r++ ) {
			EnergyMap em = pose.energies().residue_total_energies(r);

			D.push_back( em[st] );
		}

		std::string file2 = data_file_name+"._tmp_";

		write_vector_to_file(D, file2);
		//TS_ASSERT_FILE_EQ_AS_DOUBLE(data_file_name.c_str(), file2.c_str(),
		//       abs_p, rel_p);
		CxxTest::doAssertFileEQ_AsDouble(__FILE__, __LINE__,
			data_file_name.c_str(), data_file_name.c_str(),
			file2.c_str(), file2.c_str(), abs_p, rel_p, 0);

	}

	void write_vector_to_file(std::vector<double> const & v, std::string filename) {
		std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
		if ( !file ) {
			Error() << "write_vector_to_file: Unable to open file:" << filename << " for writing!!!\n";
			return;
		}
		for ( unsigned int i=0; i<v.size(); i++ ) file << v[i] << "\n";
		file.close();
	}

	void one_score_type_test_old(scoring::ScoreType st, std::string pdb_file_name, std::string utracer_file_name)
	{
		//std::cout << " Testing score: " << scoring::name_from_score_type(st) << "..." << std::endl;

		Pose pose;
		core::import_pose::pose_from_file(pose, pdb_file_name, core::import_pose::PDB_file);

		SymmetricScoreFunction scorefxn;
		scorefxn.set_weight(st, 1.0 );

		Energy score = scorefxn( pose );
		/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( pose );

		test::UTracer UT(utracer_file_name);
		UT << (int)st << " Energy=" << score << "\n";

		for ( Size r=1; r<=pose.total_residue(); r++ ) {
			EnergyMap em = pose.energies().residue_total_energies(r);

			UT << " residue: " << r << " Energy=" << em[st] << "\n";
		}
	}

	void
	test_sym_vs_asym_score()
	{
		Real score_delta_threshold = 1.2;
		Pose pose, pose_asym;
		core::import_pose::pose_from_file(pose, "core/scoring/symmetry/test_in.pdb" , core::import_pose::PDB_file);
		core::pose::symmetry::make_symmetric_pose( pose );

		// Note: hydrogens have been stripped from both pdbs so they are placed using the same logic
		core::import_pose::pose_from_file(pose_asym, "core/scoring/symmetry/symm_test.pdb" , core::import_pose::PDB_file);

		pose.dump_pdb( "test_sym_vs_asym_score_sym_pose.pdb" );
		pose_asym.dump_pdb( "test_sym_vs_asym_score_asym_pose.pdb" );


		core::scoring::ScoreFunctionOP scorefxn =
			core::scoring::ScoreFunctionFactory::create_score_function( "score13_env_hb" );

		Real sym_score = (*scorefxn)(pose);
		Real asym_score = (*scorefxn)(pose_asym);
		//std::cout << " Comparing the symmetrical score: " << sym_score << " to the assymetric: " <<
		//     asym_score << " difference: " << sym_score - asym_score << std::endl;
		TS_ASSERT( std::fabs( sym_score - asym_score ) < score_delta_threshold );
	}
};
