// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreTest.cxxtest.hh
/// @brief  unified scoring test.
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project headers
//#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.ScoreTest.cxxtest");

// using declarations
using namespace core;
using namespace scoring;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;

///////////////////////////////////////////////////////////////////////////
/// @name ScoreTest
/// @brief: unified tests for difference score functions/methods
///////////////////////////////////////////////////////////////////////////
class ScoreTest : public CxxTest::TestSuite {

public:
	//chemical::ResidueTypeSetCAP residue_set;

	void setUp() {
		core_init_with_additional_options( "-no_optH" );
		//residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	void tearDown() {}

	void test_ScoreTest() {
		//one_score_type_test(, "core/scoring/test_in.pdb", "core/scoring/.u");

		one_score_type_test(fa_atr, "core/scoring/test_in.pdb", "core/scoring/fa_atr.data");

		one_score_type_test(fa_rep, "core/scoring/test_in.pdb", "core/scoring/fa_rep.data");
		one_score_type_test(fa_sol, "core/scoring/test_in.pdb", "core/scoring/fa_sol.data");

		one_score_type_test(fa_intra_atr, "core/scoring/test_in.pdb", "core/scoring/fa_intra_atr.data");
		one_score_type_test(fa_intra_rep, "core/scoring/test_in.pdb", "core/scoring/fa_intra_rep.data");
		one_score_type_test(fa_intra_sol, "core/scoring/test_in.pdb", "core/scoring/fa_intra_sol.data");

		/*
		one_score_type_test(coarse_fa_atr, "core/scoring/test_in.pdb", "core/scoring/coarse_fa_atr.u");
		coarse_fa_rep,
		coarse_fa_sol,
		coarse_beadlj,
		*/

		one_score_type_test(mm_twist, "core/scoring/test_in.pdb", "core/scoring/mm_twist.data");
		//one_score_type_test(mm_bend, "core/scoring/test_in.pdb", "core/scoring/mm_bend.u");
		//one_score_type_test(mm_stretch, "core/scoring/test_in.pdb", "core/scoring/mm_stretch.u");
		one_score_type_test(fa_elec, "core/scoring/test_in.pdb", "core/scoring/fa_elec.data");
		one_score_type_test(atom_pair_constraint, "core/scoring/test_in.pdb", "core/scoring/atom_pair_constraint.data");
		one_score_type_test(coordinate_constraint, "core/scoring/test_in.pdb", "core/scoring/coordinate_constraint.data");

		// 0 one_score_type_test(angle_constraint, "core/scoring/test_in.pdb", "core/scoring/angle_constraint.u");
		// 0 one_score_type_test(dihedral_constraint, "core/scoring/test_in.pdb", "core/scoring/dihedral_constraint.u");
		// 0 one_score_type_test(dna_bp, "core/scoring/test_in.pdb", "core/scoring/dna_bp.u");
		// 0 one_score_type_test(dna_bs, "core/scoring/test_in.pdb", "core/scoring/dna_bs.u");

		// nw one_score_type_test(fa_pair, "core/scoring/test_in.pdb", "core/scoring/fa_pair.u");
		//one_score_type_test(fa_plane, "core/scoring/test_in.pdb", "core/scoring/fa_plane.u");
		one_score_type_test(hbond_sr_bb, "core/scoring/test_in.pdb", "core/scoring/hbond_sr_bb.data");
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
	//

	chainbreak,
	*/

	}

	void one_score_type_test(
		scoring::ScoreType st,
		std::string pdb_file_name,
		std::string data_file_name,
		double abs_p=0.0001,
		double rel_p=0.0001
	)
	{
		TR << " Testing score: " << scoring::name_from_score_type(st) << "..." << std::endl;

		Pose pose;
		core::import_pose::pose_from_pdb(pose, pdb_file_name);

		ScoreFunction scorefxn;
		scorefxn.set_weight(st, 1.0 );

		Energy score = scorefxn( pose );
		/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( pose );

		std::vector<double> D;

		D.push_back( score );

		for(Size r=1; r<=pose.total_residue(); r++) {
			EnergyMap em = pose.energies().residue_total_energies(r);

			D.push_back( em[st] );
		}

		std::string file2 = data_file_name+"._tmp_";

		write_vector_to_file(D, file2);
		//TS_ASSERT_FILE_EQ_AS_DOUBLE(data_file_name.c_str(), file2.c_str(),
		//							abs_p, rel_p);
		CxxTest::doAssertFileEQ_AsDouble(__FILE__, __LINE__,
										 data_file_name.c_str(), data_file_name.c_str(),
										 file2.c_str(), file2.c_str(), abs_p, rel_p, 0);

	}

	void write_vector_to_file(std::vector<double> const & v, std::string filename) {
		std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
		if(!file) {
			Error() << "write_vector_to_file: Unable to open file:" << filename << " for writing!!!\n";
			return;
		}
		for(unsigned int i=0; i<v.size(); i++) file << v[i] << "\n";
		file.close();
	}

	void one_score_type_test_old(
		scoring::ScoreType st,
		std::string pdb_file_name,
		std::string utracer_file_name
	)
	{
		TR << " Testing score: " << scoring::name_from_score_type(st) << "..." << std::endl;

		Pose pose;
		core::import_pose::pose_from_pdb(pose, pdb_file_name);

		ScoreFunction scorefxn;
		scorefxn.set_weight(st, 1.0 );

		Energy score = scorefxn( pose );
		/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( pose );

		test::UTracer UT(utracer_file_name);
		UT << (int)st << " Energy=" << score << "\n";

		for(Size r=1; r<=pose.total_residue(); r++) {
			EnergyMap em = pose.energies().residue_total_energies(r);

			UT << " residue: " << r << " Energy=" << em[st] << "\n";
		}
	}

};
