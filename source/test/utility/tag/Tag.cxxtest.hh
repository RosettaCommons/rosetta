// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/tag/lexically_parse_tag.cxxtest.hh
/// @brief test lexically parsing simple xml-like files
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <sstream>
#include <string>


class TagTests : public CxxTest::TestSuite {

public:

	void setUp() {
	}

	void do_test_one_tag_fail(
		std::string const & in_tag,
		platform::Size col
	) {
		//std::cout << in_tag << std::endl;
		std::stringstream in;
		in << in_tag;
		try {
			utility::tag::TagCOP tags(utility::tag::Tag::create(in));
			TS_ASSERT(false);
		} catch ( utility::excn::EXCN_BadInput e ){
			std::stringstream expected_error;
			if(col==0){
				expected_error
					<< "Tag::read - parse error, printing backtrace." << std::endl
					<< std::endl;
			} else {
				expected_error
					<<"Tag::read - parse error, printing backtrace." << std::endl
					<< std::endl
					<< "Tag::read - parse error - file:istream line:1 column:" << col << " - " << in_tag << std::endl
					<< "Tag::read - parse error - file:istream line:1 column:" << col << " -" << std::string(col, ' ') << "^" << std::endl
					<< std::endl;
			}
			if( e.msg() != expected_error.str() ) {
				std::cout << "expected error: '" << expected_error.str() << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error.str() );
		}
	}

	void test_fail_1() {
		do_test_one_tag_fail("<>", 1);
	}

	void test_fail_2() {
		do_test_one_tag_fail("</>", 1);
	}

	void test_fail_3() {
		do_test_one_tag_fail("<abc=def/>", 1);
	}

	void test_fail_4() {
		do_test_one_tag_fail("<abc>", 1);
	}

	void test_fail_5() {
		do_test_one_tag_fail("<", 1);
	}

	void test_fail_6() {
		do_test_one_tag_fail(">", 1);
	}

	void test_fail_7() {
		do_test_one_tag_fail("><", 1);
	}

	void test_fail_8() {
		do_test_one_tag_fail("<<", 1);
	}

	void test_fail_9() {
		do_test_one_tag_fail(">>", 1);
	}

	void test_fail_10() {
		do_test_one_tag_fail("<?>", 1);
	}

	void test_fail_12() {
		do_test_one_tag_fail("<?xml>", 1);
	}

	void test_55() {
		// this should fail because there is no actual tag that is defined
		do_test_one_tag_fail("<?xml?>", 1);
	}

	void test_fail_13() {
		do_test_one_tag_fail("<!--", 1);
	}

	void test_fail_14() {
		do_test_one_tag_fail("<!-- -- -->", 1);
	}

	void test_fail_15() {
		do_test_one_tag_fail("<!-->", 1);
	}

	void test_56() {
		// this should fail because the <?xml ?> tag should come first
		do_test_one_tag_fail("<name abc=def/> <?xml ?> ",0);
	}

	void test_57() {
		// this should fail because the <?xml ?> tag should come first
		do_test_one_tag_fail("<name abc=def/> <?xml bla bla bla ?>",0);
	}

	void test_58() {
		// this should fail because the <?xml ?> tag should come first
		do_test_one_tag_fail("<name abc=def/> <?xml <name abc=hij> ?>",0);
	}

	void test_59() {
		do_test_one_tag_fail("<name abc=def/> <?xml <!-- ---> ?>",0);
	}


	void do_test_one_tag(
		std::string const & in_tag,
		std::string const & key,
		std::string const & val) {
		//std::cout << in_tag << std::endl;
		std::stringstream in;
		in << in_tag;
		utility::tag::TagCOP tags(utility::tag::Tag::create(in));
		TS_ASSERT_EQUALS(tags->getOption<std::string>(key), val);
		TS_ASSERT_EQUALS(tags->size(), 1);
	}

	void test_4() {
		do_test_one_tag("<name abc=def/>", "abc", "def");
	}

	void test_5() {
		do_test_one_tag("abc def <name abc=def/>", "abc", "def");
	}

	void test_6() {
		do_test_one_tag("<name abc=\" \"/>", "abc", " ");
	}

	void test_7() {
		do_test_one_tag(">>> <name abc=\" \"/> >>>", "abc", " ");
	}

	void test_8() {
		do_test_one_tag("<?xml ?> data etc. <name abc=\" \"/> >>>", "abc", " ");
	}

	void test_17() {
		do_test_one_tag("<name abc=\" <?xml \"/>", "abc", " <?xml ");
	}

	void test_18() {
		do_test_one_tag("<name abc=\" <!-- \"/>", "abc", " <!-- ");
	}

	void test_19() {
		do_test_one_tag("<name abc=\" < \"/>", "abc", " < ");
	}

	void test_20() {
		do_test_one_tag("<name abc=\" ?> \"/>", "abc", " ?> ");
	}

	void test_21() {
		do_test_one_tag("<name abc=\" --> \"/>", "abc", " --> ");
	}

	void test_22() {
		do_test_one_tag("<name abc=\" ---> \"/>", "abc", " ---> ");
	}

	void test_23() {
		do_test_one_tag("<name abc=\" -- \"/>", "abc", " -- ");
	}

	void test_24() {
		do_test_one_tag("<name abc=\" > \"/>", "abc", " > ");
	}

	void test_25() {
		do_test_one_tag("<name abc=\" ' \"/>", "abc", " ' ");
	}

	void test_26() {
		do_test_one_tag("<name abc=\" '' \"/>", "abc", " '' ");
	}

	void test_40() {
		do_test_one_tag("<name abc=\" <?xml \n\"/>", "abc", " <?xml \n");
	}

	void test_41() {
		do_test_one_tag("<name abc=\" <!-- \n\"/>", "abc", " <!-- \n");
	}

	void test_42() {
		do_test_one_tag("<name abc=\" < \n\"/>", "abc", " < \n");
	}

	void test_43() {
		do_test_one_tag("<name abc=\" ?> \n\"/>", "abc", " ?> \n");
	}

	void test_44() {
		do_test_one_tag("<name abc=\" --> \n\"/>", "abc", " --> \n");
	}

	void test_45() {
		do_test_one_tag("<name abc=\" ---> \n\"/>", "abc", " ---> \n");
	}

	void test_46() {
		do_test_one_tag("<name abc=\" -- \n\"/>", "abc", " -- \n");
	}

	void test_47() {
		do_test_one_tag("<name abc=\" > \n\"/>", "abc", " > \n");
	}

	void test_48() {
		do_test_one_tag("<name abc=\" ' \n\"/>", "abc", " ' \n");
	}

	void test_49() {
		do_test_one_tag("<name abc=\" '' \n\"/>", "abc", " '' \n");
	}

	void test_66() {
		do_test_one_tag("<name abc=def/> <!-- -->", "abc", "def");
	}

	void test_60() {
		do_test_one_tag("<?xml?> <name abc=def/>", "abc", "def");
	}

	void test_61() {
		do_test_one_tag("<?xml ?> <name abc=def/>", "abc", "def");
	}

	void test_62() {
		do_test_one_tag("<?xml bla bla bla ?> <name abc=def/>", "abc", "def");
	}

	void test_63() {
		do_test_one_tag("<?xml <name abc=hij> ?> <name abc=def/>", "abc", "def");
	}

	void test_64() {
		do_test_one_tag("<?xml <!-- --> ?> <name abc=def/>", "abc", "def");
	}

	void test_65() {
		do_test_one_tag("<name abc=def/> <!-- - -->", "abc", "def");
	}

	void test_67() {
		do_test_one_tag("<name abc=def/> <!-- bla bla bla -->", "abc", "def");
	}

	void test_68() {
		do_test_one_tag("<name abc=def/> <!-- <name abc=hij> -->", "abc", "def");
	}

	void test_69() {
		do_test_one_tag("<name abc=def/> <!-- <?xml ?> -->", "abc", "def");
	}

	void test_70() {
		do_test_one_tag(" <!-- --> <name abc=def/>", "abc", "def");
	}

	void test_71() {
		do_test_one_tag(" <!-- - --> <name abc=def/>", "abc", "def");
	}

	void test_72() {
		do_test_one_tag(" <!-- bla bla bla --> <name abc=def/>", "abc", "def");
	}

	void test_73() {
		do_test_one_tag(" <!-- <name abc=hij> --> <name abc=def/>", "abc", "def");
	}

	void test_74() {
		do_test_one_tag(" <!-- <?xml ?> --> <name abc=def/>", "abc", "def");
	}

	void test_75() {
		do_test_one_tag("<name abc=def></name>", "abc", "def");
	}

	void test_76() {
		do_test_one_tag("<name abc=def>123</name>", "abc", "def");
	}

	void test_77() {
		do_test_one_tag(" <name abc=def> </name> ", "abc", "def");
	}

	void test_78() {
		do_test_one_tag("< name abc=def > < / name > ", "abc", "def");
	}

	void test_79() {
		do_test_one_tag("<name abc=def> <!-- --> </name> ", "abc", "def");
	}

	void test_80() {
		do_test_one_tag("<name abc=def> <!-- < --> </name> ", "abc", "def");
	}

	void test_81() {
		do_test_one_tag("<name abc=def> <!-- > --> </name> ", "abc", "def");
	}

	void test_82() {
		do_test_one_tag("<name abc=def> <!-- <?xml ?> --> </name> ", "abc", "def");
	}

	void test_83() {
		do_test_one_tag("<name abc=def> <!--  --> <!-- --> </name> ", "abc", "def");
	}

	void test_84() {
		do_test_one_tag("<name abc=def> >>> </name> ", "abc", "def");
	}

	void do_test_two_kv_tag(
		std::string const & in_tag,
		std::string const & key1,
		std::string const & val1,
		std::string const & key2,
		std::string const & val2) {
		//std::cout << in_tag << std::endl;
		std::stringstream in;
		in << in_tag;
		utility::tag::TagCOP tags(utility::tag::Tag::create(in));
		TS_ASSERT_EQUALS(tags->getOption<std::string>(key1), val1);
		TS_ASSERT_EQUALS(tags->getOption<std::string>(key2), val2);
		TS_ASSERT_EQUALS(tags->size(), 1);
	}

	void test_85() {
		do_test_two_kv_tag("<name abc=def xyz=uvw/>", "abc", "def", "xyz", "uvw");
	}

	void test_86() {
		do_test_two_kv_tag("< name abc = def xyz = uvw / >", "abc", "def", "xyz", "uvw");
	}

	void test_87() {
		do_test_two_kv_tag("< name abc =\"def\" xyz = \"uvw\"/>", "abc", "def", "xyz", "uvw");
	}


	void do_test_nested_tags(
		std::string const & in_tag,
		std::string const & key1,
		std::string const & val1,
		std::string const & key2,
		std::string const & val2) {
		//std::cout << in_tag << std::endl;
		std::stringstream in;
		in << in_tag;
		utility::tag::TagCOP tags(utility::tag::Tag::create(in));
		TS_ASSERT_EQUALS(tags->getOption<std::string>(key1), val1);
		TS_ASSERT_EQUALS(tags->size(), 2);
		for(
			utility::tag::Tag::tags_t::const_iterator
				tag_iter = tags->getTags().begin(),
				tag_iter_end = tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter
		) {
			TS_ASSERT_EQUALS((*tag_iter)->getOption<std::string>(key2), val2);
			TS_ASSERT_EQUALS((*tag_iter)->size(), 1);
		}
	}

	void test_88() {
		do_test_nested_tags("<name abc=def> <name2 xyz=uvw/> </name>", "abc", "def", "xyz", "uvw");
	}

	void test_89() {
		do_test_nested_tags("<name abc=def>123<name2 xyz=uvw/>456</name>", "abc", "def", "xyz", "uvw");
	}

	void test_90() {
		do_test_nested_tags("<name abc=def> >> <name2 xyz=uvw/> <!-- <name3 hij=klm/> --> </name>", "abc", "def", "xyz", "uvw");
	}


};
