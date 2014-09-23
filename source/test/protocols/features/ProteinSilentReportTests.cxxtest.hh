// -*- mode:c++;tab-width:2;indent-tabs-mode:nil;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinSilentReportTests.cxxtest.hh
/// @brief  Test Protein Silent Report
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/features/ProteinSilentReport.hh>
#include <protocols/features/util.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/file/file_sys_util.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("protocols.features.ProteinSilentReportTests.cxxtest");

class ProteinSilentReportTests : public CxxTest::TestSuite {

public:

	void
	setUp() {
    using protocols::features::ProteinSilentReport;
    using core::pose::tag_into_pose;
    using utility::sql_database::DatabaseSessionManager;

    std::string db_fname("protein_silent_report.db3");

		core_init();
		pose_1ten_ = fullatom_pose_from_string( pdb_string_1ten() );
    tag_into_pose( pose_1ten_, "1ten");

    core::scoring::ScoreFunctionOP scfxn(core::scoring::get_score_function());
    scfxn->score(pose_1ten_);

    protein_silent_report_ = protocols::features::ProteinSilentReportOP( new ProteinSilentReport() );

    utility::file::file_delete(db_fname);
		db_session_ = basic::database::get_db_session(db_fname);
	}

	void test_main() {
		write_schema();
    write_pose_to_db();
    read_pose_from_db();
	}

	void write_schema() {
    protein_silent_report_->write_schema_to_db(db_session_);
	}

	void write_pose_to_db() {
    protein_silent_report_->apply(pose_1ten_, db_session_);
	}

  void read_pose_from_db() {
    using core::scoring::CA_rmsd;
    using core::scoring::all_atom_rmsd;
    using core::scoring::CA_gdtmm;
    using protocols::features::StructureID;

		TR << "retrieving struct_ids from DB" << std::endl;
		utility::vector1<StructureID> struct_ids = protocols::features::struct_ids_from_tag(db_session_,"1ten");
		TS_ASSERT(struct_ids.size() == 1);

    core::pose::Pose copy_1ten;
    protein_silent_report_->load_pose(db_session_, struct_ids[1], copy_1ten);

    TS_ASSERT(pose_1ten_.total_residue() == copy_1ten.total_residue());
    TS_ASSERT(pose_1ten_.sequence() == copy_1ten.sequence());
    TS_ASSERT(pose_1ten_.annotated_sequence() == copy_1ten.annotated_sequence());
    TS_ASSERT(pose_1ten_.is_fullatom() == copy_1ten.is_fullatom());
    TS_ASSERT(pose_1ten_.fold_tree() == copy_1ten.fold_tree());
    TS_ASSERT(pose_1ten_.conformation().chain_endings() == copy_1ten.conformation().chain_endings());

    for(core::Size i = 1; i <= pose_1ten_.total_residue(); ++i){
      core::conformation::Residue const & orig = pose_1ten_.residue(i);
      core::conformation::Residue const & copy = copy_1ten.residue(i);

      TS_ASSERT(pose_1ten_.secstruct(i) == copy_1ten.secstruct(1));
      TS_ASSERT_DELTA(orig.mainchain_torsion(1), copy.mainchain_torsion(1),0.000001);
      TS_ASSERT_DELTA(orig.mainchain_torsion(2), copy.mainchain_torsion(2),0.000001);
      TS_ASSERT_DELTA(orig.mainchain_torsion(3), copy.mainchain_torsion(3),0.000001);
      if (orig.nchi() >= 1) TS_ASSERT_DELTA(orig.chi(1), copy.chi(1),0.000001);
      if (orig.nchi() >= 2) TS_ASSERT_DELTA(orig.chi(2), copy.chi(2),0.000001);
      if (orig.nchi() >= 3) TS_ASSERT_DELTA(orig.chi(3), copy.chi(3),0.000001);
      if (orig.nchi() >= 4) TS_ASSERT_DELTA(orig.chi(4), copy.chi(4),0.000001);


      TS_ASSERT(orig.natoms() ==  copy.natoms());
      for(core::Size atom = 1; atom <= orig.natoms(); ++atom)
      {
    	  core::Vector orig_coord(orig.xyz(atom));
    	  core::Vector copy_coord(copy.xyz(atom));
    	  TS_ASSERT_DELTA(orig_coord.x(),copy_coord.x(),0.000001);
    	  TS_ASSERT_DELTA(orig_coord.y(),copy_coord.y(),0.000001);
    	  TS_ASSERT_DELTA(orig_coord.z(),copy_coord.z(),0.000001);

      }

    }

    TS_ASSERT_DELTA(CA_rmsd(pose_1ten_, copy_1ten),0.0,0.000001);
    TS_ASSERT_DELTA(all_atom_rmsd(pose_1ten_,copy_1ten),0.0,0.000001);
    TS_ASSERT_DELTA(CA_gdtmm(pose_1ten_, copy_1ten),1.0,0.000001);


  }

private:
	core::pose::Pose pose_1ten_;
  protocols::features::ProteinSilentReportOP protein_silent_report_;
  utility::sql_database::sessionOP db_session_;
};
