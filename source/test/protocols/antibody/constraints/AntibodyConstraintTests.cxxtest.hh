// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/constraints/AntibodyConstraintTests.cxxtest.hh
/// @brief  tests for the AntibodyInfo class
/// @author Jared Adolf-Bryfogle


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMover.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <utility/vector1.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

#define BFE BOOST_FOREACH
using namespace protocols::antibody;
using namespace protocols::antibody::constraints;

using utility::vector1;

static thread_local basic::Tracer TR("protocols.antibody.AntibodyConstraintTests");


class AntibodyConstraintTests : public CxxTest::TestSuite {
	core::pose::Pose pose;
	AntibodyInfoOP ab_info;
	ParatopeSiteConstraintMoverOP paratope_mover;
	ParatopeEpitopeSiteConstraintMoverOP para_epitope_mover;
	CDRDihedralConstraintMoverOP cdr_dihedral_cst_mover;


public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_pdb(pose, "protocols/antibody/1bln_AB_aho.pdb"); //AHO renumbered pose
		ab_info = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );
		paratope_mover = ParatopeSiteConstraintMoverOP( new constraints::ParatopeSiteConstraintMover(ab_info) );
		para_epitope_mover = ParatopeEpitopeSiteConstraintMoverOP( new constraints::ParatopeEpitopeSiteConstraintMover(ab_info) );
		cdr_dihedral_cst_mover = CDRDihedralConstraintMoverOP( new CDRDihedralConstraintMover(ab_info));

	}

	void test_constraints(){
		cdr_dihedral_cst_mover->set_cdr(l1);
		TS_ASSERT_THROWS_NOTHING(cdr_dihedral_cst_mover->apply(pose));

		core::scoring::constraints::ConstraintSetCOP  csts = pose.constraint_set();
		TS_ASSERT(csts->has_constraints());
		TS_ASSERT(protocols::antibody::constraints::cdr_has_res_constraints(ab_info, pose, l1, "Dihedral"))
			utility::vector1<core::scoring::constraints::ConstraintCOP> cst_list = csts->get_all_constraints();
		TS_ASSERT(cst_list.size() == 32); //16 residues, one constraint for each dihedral on the CDR.
		TS_ASSERT(cst_list[1]->type() == "Dihedral");
		pose.remove_constraints();

		TS_ASSERT_THROWS_NOTHING(protocols::antibody::constraints::add_harmonic_dihedral_cst_general(ab_info, pose, l1, 23.0, 43.0));

		csts = pose.constraint_set();
		TS_ASSERT(csts->has_constraints());
		cst_list = csts->get_all_constraints();
		TS_ASSERT(cst_list.size() == 32); //16 residues, one constraint for each dihedral on the CDR.
		TS_ASSERT(cst_list[1]->type() == "Dihedral");
		pose.remove_constraints();

		//protocols::antibody::constraints::ParatopeEpitopeSiteConstraintMoverOP
	}
};
