// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/tcr/TCRmodelSeqParserTests
/// @brief  tests for TCR seqence parsing, regular expresssion, cdr assignment
/// @author Ragul Gowthaman

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
// Core Headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
// Protocol Headers
#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/grafting_util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.tcr.TCRmodelGraftingTests");
class TCRmodelGraftingTests : public CxxTest::TestSuite {

	protocols::tcr::TCRseqInfo::tcrposi a_aho_posi;
	protocols::tcr::TCRseqInfo::tcrposi b_aho_posi;
	protocols::tcr::TCRmodel::tcrtmplts tcrseg;
	std::string tmplt_pdb_path;

public:
	void setUp(){
		core_init();
		tcrseg.fr.tid = "4x6b_A";
		tcrseg.cdr3.tid = "4x6b_A";
		tmplt_pdb_path =  "protocols/tcr/4x6b_aho.pdb";
	}

	void tearDown(){
		tcrseg.fr.tpiece->clear();
		tcrseg.cdr3.tpiece->clear();
	}

	void test_tcr_grafting(){
		using namespace protocols::tcr;
		initialize_aho_numbers(a_aho_posi, b_aho_posi);
		graft_framework(tcrseg.fr, tmplt_pdb_path, a_aho_posi);
		graft_cdr(tcrseg.cdr3, tmplt_pdb_path, a_aho_posi.cdr3.begin, a_aho_posi.cdr3.end);

		//expected
		std::string fr_tmplt_seq = "VEQDPGPLSVPEGAIVSLNCTYSNSAFQYFMWYRQYSRKGPELLMYTYSSGNKEDGRFTAQVDKSSKYISLFIRDSQPSDSATYLCAMSTSLPNAGKSTFGDGTTLTVK";
		std::string cdr3a_tmplt_seq = "AMSTSLPNAGKST";

		TR << "Testing TCRmodelGraftingTests" << std::endl;
		TR << "Grafted frwmework sequence : " << fr_tmplt_seq << " " << tcrseg.fr.tpiece->sequence() <<std::endl;
		TR << "Grafted cdr3 sequence : " << cdr3a_tmplt_seq << " " << tcrseg.cdr3.tpiece->sequence() <<std::endl;
		TS_ASSERT_EQUALS(fr_tmplt_seq, tcrseg.fr.tpiece->sequence());
		TS_ASSERT_EQUALS(cdr3a_tmplt_seq, tcrseg.cdr3.tpiece->sequence());
	}

};


