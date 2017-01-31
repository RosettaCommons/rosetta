// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/hybridization/TMalign.cxxtest.hh
/// @brief  test suite for TMalign functions in hybridization
/// @details Unit test sets up 2 poses, perturbs one, then uses TMalign superimpose. This should yield a TMscore of 1, although the test is only set to .95.
/// @author Steve Bertolani sjbertolani@ucdavis.edu

// Test headers
#include <cxxtest/TestSuite.h>

// c++ headers
#include <cmath>

// initialization headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

//core
#include <core/pose/Pose.hh>
#include <core/pose/util.hh> // for initialize_atommap

//basic
#include <basic/Tracer.hh>

//utility
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// TMalign headers
#include <protocols/hybridization/TMalign.hh>
#include <protocols/hybridization/util.hh>


static basic::Tracer TR("test.protocols.hybridization.TMalign");

class hybridization_TMalign_Tests : public CxxTest::TestSuite {
public:
	hybridization_TMalign_Tests() {}

	void setUp() {

		core_init();

		/* Setup a query pose and reference pose, perturb the query pose, then the test function will re-superimpose them (they are from the same PDB lines, so they had better line up again */


		TR << "Running setup for TMalign unit test" << std::endl;;
		query_pose_ = create_test_in_pdb_pose(); // 116 residues
		TR << "Created query pose" << std::endl;
		start_ = query_pose_.residue(1).xyz("CA");  //store a xyz vector for res1 to make sure the perturbation works

		ref_pose = create_test_in_pdb_pose();
		TR << "Created reference pose" << std::endl;

		TR << "Perturbing the query pose" << std::endl;

		for ( core::Size i=1; i <= query_pose_.size(); ++i ) {
			for ( core::Size j=1; j <=query_pose_.residue(i).natoms(); ++j ) {
				query_pose_.set_xyz(core::id::AtomID(j,i), query_pose_.residue(i).atom(j).xyz()+core::Vector(1,1,1));
			}
		}
		perturb_ = query_pose_.residue(1).xyz("CA");  // perturbed xyz vector to compare

	}

	void tearDown() {}

	void test_align_pose1_onto_pose2() {


		std::list <core::Size> pose_residue_list; // all residue numbers in ipose, used for transformation after alignment
		for ( core::Size ires = 1; ires <= query_pose_.size(); ++ires ) {
			pose_residue_list.push_back(ires);
		}

		std::list <core::Size> ref_pose_residue_list; // all residue numbers in ipose, used for transformation after alignment
		for ( core::Size ires = 1; ires <= ref_pose.size(); ++ires ) {
			ref_pose_residue_list.push_back(ires);
		}

		protocols::hybridization::TMalign tm_align;
		tm_align.apply(query_pose_,ref_pose,pose_residue_list,ref_pose_residue_list);
		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, query_pose_, core::id::BOGUS_ATOM_ID );
		core::Size n_mapped_residues=0;

		tm_align.alignment2AtomMap(query_pose_, ref_pose, pose_residue_list,ref_pose_residue_list, n_mapped_residues, atom_map);
		core::Size normalize_length = query_pose_.size() < ref_pose.size() ? query_pose_.size() : ref_pose.size();
		core::Real TMscore = tm_align.TMscore(normalize_length);
		utility::vector1< core::Real > aln_cutoffs;
		aln_cutoffs.push_back(2);
		aln_cutoffs.push_back(1.5);
		aln_cutoffs.push_back(1.0);
		aln_cutoffs.push_back(0.5);
		core::Real min_coverage = 0.2;
		protocols::hybridization::partial_align(query_pose_,ref_pose, atom_map, pose_residue_list, true, aln_cutoffs, min_coverage);
		TR << "TMaligned with a TMscore of " << TMscore << std::endl;
		end_ = query_pose_.residue(1).xyz("CA");  // should approx recover so that the difference very small between start_ and end_

		TR << "start_ " << start_.to_string() << std::endl;
		TR << "perturb_ " << perturb_.to_string() << std::endl;
		TR << "end_ " << end_.to_string() << std::endl;

		TS_ASSERT( TMscore >= 0.95 );
		TS_ASSERT( start_ - perturb_ != core::Vector(0,0,0) );  //check perturbation moved the protein
		TS_ASSERT_DELTA( start_.normalized().dot(end_.normalized()) /* equals 1 when it works */, 1, .05 ); // normalized will equal 1, but just to be sure just guarantee the delta is small


		// check the getters for the kabsch t and u
		//  core::Vector(-1,-1,-1) // this is what we perturbed by above
		numeric::xyzVector< core::Real> testvector(-1,-1,-1);
		numeric::xyzVector< core::Real> testvector2(-1,-1,-1);
		TR.Debug << "testvector" << std::endl;
		TR.Debug << testvector.to_string() << std::endl;

		// should be -1,-1,-1 .. but not actually (use debug and gdb>p tm_align.get_t() )
		// compare with gdb>p testvector
		TR.Debug << " T vector" << tm_align.get_t().to_string() << std::endl;
		TR.Debug << testvector.dot( tm_align.get_t() ) << std::endl;;

		TS_ASSERT( std::abs( testvector.dot( tm_align.get_t() ) - 3 ) < 0.05);

		TR.Debug << " U matrix" << std::endl;
		tm_align.get_u().show(TR); // << std::endl;

		//

	}
private:
	core::pose::Pose query_pose_;
	core::pose::Pose ref_pose;
	core::Vector start_;
	core::Vector perturb_;
	core::Vector end_;
}; // class hybridization_TMalign_Tests
