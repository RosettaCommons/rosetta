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
// Protocol Headers
#include <protocols/antibody/grafting/util.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/util.hh>
#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.tcr.TCRmodelSeqParserTests");

class TCRmodelSeqParserTests : public CxxTest::TestSuite {

	std::string aseq;
	std::string bseq;
	std::string cdr1a_num;
	std::string cdr2a_num;
	std::string cdr3a_num;
	std::string cdr1b_num;
	std::string cdr2b_num;
	std::string cdr3b_num;
	protocols::tcr::TCRseqInfo::tcrsegs atcr_re;
	protocols::tcr::TCRseqInfo::tcrsegs btcr_re;
	protocols::tcr::TCRseqInfo::tcrsegs atcr_num;
	protocols::tcr::TCRseqInfo::tcrsegs btcr_num;
	protocols::tcr::TCRseqInfo::tcrposi a_aho_posi;
	protocols::tcr::TCRseqInfo::tcrposi b_aho_posi;
	protocols::tcr::TCRseqInfo::tcrposi atcr_posi;
	protocols::tcr::TCRseqInfo::tcrposi btcr_posi;
	protocols::tcr::TCRseqInfo::tcrposi anum_posi;
	protocols::tcr::TCRseqInfo::tcrposi bnum_posi;

public:

	void setUp(){
		core_init();
		aseq = "QSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGLQLLLKYYSGDPVVQGVNGFEAEFSKSNSSFHLRKASVHWSDSAVYFCAVSGFASALTFGSGTKVIVLPYIQNPEPAVYALKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDATVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPC";
		bseq = "EAAVTQSPRNKVAVTGGKVTLSCNQTNNHNNMYWYRQDTGHGLRLIHYSYGAGSTEKGDIPDGYKASRPSQENFSLILELATPSQTSVYFCASGGGGTLYFGAGTRLSVLEDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADC";
		cdr1a_num = "23:33";
		cdr2a_num = "47:74";
		cdr3a_num = "91:100";
		cdr1b_num = "24:33";
		cdr2b_num = "47:75";
		cdr3b_num = "92:100";
	}

	void tearDown(){
	}

	void test_tcr_seq_parsing(){

#ifdef __ANTIBODY_GRAFTING__
		//assignment by regex not compiling on some test servers (linux.clang, linux.gcc)
		//probably due to lower versions of clang/gcc
		if( ! protocols::antibody::grafting::antibody_grafting_usable() ) {
			TR << "SKIPPING TCRmodelSeqParserTests : Compiler does not have full support for C++11 regex." <<std::endl;
			return; // Don't attempt to run test on system which doesn't support it.
		}

		//expected
		std::string cdr1a_seq = "KYSYSATPYLF";
		std::string cdr3a_seq = "AVSGFASALT";
		std::string cdr1b_seq = "NQTNNHNNMY";
		std::string cdr3b_seq = "ASGGGGTLY";

		using namespace protocols::tcr;
		initialize_aho_numbers(a_aho_posi, b_aho_posi);

		//Assign CDR's using regex
		assign_achain_CDRs_using_REGEX(aseq, atcr_re, atcr_posi);
		assign_bchain_CDRs_using_REGEX(bseq, btcr_re, btcr_posi);
		adjust_position_for_chain(btcr_posi, atcr_re.truncdomain.length() );
		TR << "Testing TCRmodelSeqParserTests for regex assignment" << std::endl;
		TR << "RegEx assignment, CDR1 Alpha : " << atcr_re.cdr1 << " " << cdr1a_seq <<std::endl;
		TR << "RegEx assignment, CDR3 Alpha : " << atcr_re.cdr3 << " " << cdr3a_seq <<std::endl;
		TR << "RegEx assignment, CDR1 Beta : " << btcr_re.cdr1 << " " << cdr1b_seq <<std::endl;
		TR << "RegEx assignment, CDR3 Beta : " << btcr_re.cdr3 << " " << cdr3b_seq <<std::endl;
		TS_ASSERT_EQUALS(atcr_re.cdr1, cdr1a_seq);
		TS_ASSERT_EQUALS(atcr_re.cdr3, cdr3a_seq);
		TS_ASSERT_EQUALS(btcr_re.cdr1, cdr1b_seq);
		TS_ASSERT_EQUALS(btcr_re.cdr3, cdr3b_seq);

		//Assign CDR's using user numbering
		//start/end positions provided by user
		anum_posi.cdr1 = string_to_CDRbounds(cdr1a_num);
		anum_posi.cdr2hv4 = string_to_CDRbounds(cdr2a_num);
		anum_posi.cdr3 = string_to_CDRbounds(cdr3a_num);
		bnum_posi.cdr1 = string_to_CDRbounds(cdr1b_num);
		bnum_posi.cdr2hv4 = string_to_CDRbounds(cdr2b_num);
		bnum_posi.cdr3 = string_to_CDRbounds(cdr3b_num);
		assign_CDRs_using_numbers(aseq, anum_posi.cdr1, anum_posi.cdr2hv4, anum_posi.cdr3, a_aho_posi.cap, atcr_num, atcr_posi);
		assign_CDRs_using_numbers(bseq, bnum_posi.cdr1, bnum_posi.cdr2hv4, bnum_posi.cdr3, b_aho_posi.cap, btcr_num, btcr_posi);
		adjust_position_for_chain(btcr_posi, atcr_num.truncdomain.length() );
		TR << "Testing TCRmodelSeqParserTests for assignment by user numbering" << std::endl;
		TR << "Num assignment, CDR1 Alpha : " << atcr_num.cdr1 << " " << cdr1a_seq <<std::endl;
		TR << "Num assignment, CDR3 Alpha : " << atcr_num.cdr3 << " " << cdr3a_seq <<std::endl;
		TR << "Num assignment, CDR1 Beta : " << btcr_num.cdr1 << " " << cdr1b_seq <<std::endl;
		TR << "Num assignment, CDR3 Beta : " << btcr_num.cdr3 << " " << cdr3b_seq <<std::endl;
		TS_ASSERT_EQUALS(atcr_num.cdr1, cdr1a_seq);
		TS_ASSERT_EQUALS(atcr_num.cdr3, cdr3a_seq);
		TS_ASSERT_EQUALS(btcr_num.cdr1, cdr1b_seq);
		TS_ASSERT_EQUALS(btcr_num.cdr3, cdr3b_seq);

#endif // __ANTIBODY_GRAFTING__
	}


};


