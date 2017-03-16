// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/FragmentCandidatesTests.cxxtest.hh
/// @brief
/// @author Domini Gront

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/DiversifyCrmsdSelector.hh>
#include <protocols/frag_picker/DiversifyDihedralsSelector.hh>
#include <protocols/frag_picker/DiversifyCrmsdByClustering.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

//#include <utility/stream_util.hh>

#include <iostream>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer tr("protocols.frag_picker.FragmentCandidatesTests.cxxtest");

using namespace core;
using namespace protocols::frag_picker;

class FragmentCandidatesTests : public CxxTest::TestSuite {
public:

	FragmentCandidatesTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// test for reading fragments and vall files
	void test_read_vall_fragments() {

		chunks_ = protocols::frag_picker::VallProviderOP( new VallProvider() );
		chunks_->vallChunksFromLibrary("protocols/frag_picker/1a32A-H1.vall");
		TS_ASSERT( chunks_->size() == 229);

		in_frags_ =
			read_fragment_candidates("protocols/frag_picker/1a32A-H1.9mers",chunks_);
		TS_ASSERT( in_frags_.size() == 250);
	}

	// test for reading vall file twice (should be appended)
	void test_read_vall2() {

		VallProviderOP chunks( new VallProvider() );
		chunks->vallChunksFromLibrary("protocols/frag_picker/1a32A-H1.vall");
		TS_ASSERT( chunks->size() == 229);
		Real k1 = (Real) chunks->at(chunks->size())->key();
		tr<<"Key of the last chunk: "<<chunks->at(chunks->size())->key()<<std::endl;
		chunks->vallChunksFromLibrary("protocols/frag_picker/1a32A-H1.vall");
		TS_ASSERT( chunks->size() == 458);
		tr<<"Key of the last chunk: "<<chunks->at(chunks->size())->key()<<std::endl;
		Real k2 = (Real) chunks->at(chunks->size())->key();
		TS_ASSERT(k2 / k1 > 0.98);
		TS_ASSERT(k2 / k1 < 2.02);
			}

			void test_DiversifyCrmsdSelector() {

			utility::vector1<std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> > fs_in;
		utility::vector1<std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> > fs_out;
		scores::FragmentScoreMapOP empty_score_map( new scores::FragmentScoreMap(0) );
		for ( Size i=1; i<=in_frags_.size(); i++ ) {
			fs_in.push_back( std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP>(in_frags_[i],empty_score_map) );
		}

		FragmentSelectingRuleOP sel( new DiversifyCrmsdSelector(10,0.8) );
		sel->select_fragments( fs_in, fs_out);
		TS_ASSERT( fs_out.size() <=10);
	}

	void test_DiversifyDihedralsSelector() {

		utility::vector1<std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> > fs_in;
		utility::vector1<std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> > fs_out;
		scores::FragmentScoreMapOP empty_score_map( new scores::FragmentScoreMap(0) );
		for ( Size i=1; i<=in_frags_.size(); i++ ) {
			fs_in.push_back( std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP>(in_frags_[i],empty_score_map) );
		}

		DiversifyDihedralsSelector* sel = new DiversifyDihedralsSelector(10,20.0);
		sel->select_fragments( fs_in, fs_out);
		TS_ASSERT( fs_out.size() <=10);

		Real min_rms = 10000000.0;
		for ( Size i=2; i<=fs_out.size(); i++ ) {
			for ( Size j=1; j<i; j++ ) {
				Real r = sel->dihedral_rmsd(fs_out[i].first, fs_out[j].first);
				if ( min_rms>r ) min_rms = r;
			}
		}
		tr<<"Min dihedral rms = "<<min_rms<<" based on "<<" fs_out.size() selected fragments"<<std::endl;
		TS_ASSERT( min_rms > 20.0);
		delete sel;
	}

	void test_DiversifyCrmsdByClustering() {

		utility::vector1<std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> > fs_in;
		utility::vector1<std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> > fs_out;
		scores::FragmentScoreMapOP empty_score_map( new scores::FragmentScoreMap(0) );
		for ( Size i=1; i<=in_frags_.size(); i++ ) {
			fs_in.push_back( std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP>(in_frags_[i],empty_score_map) );
		}

		FragmentSelectingRuleOP sel( new DiversifyCrmsdByClustering(10) );
		sel->select_fragments( fs_in, fs_out);
		TS_ASSERT( fs_out.size() <=10);
	}


	// Shared finalization goes here.
	void tearDown() {}

private:
	protocols::frag_picker::VallProviderOP chunks_;
	utility::vector1<FragmentCandidateOP> in_frags_;

};
