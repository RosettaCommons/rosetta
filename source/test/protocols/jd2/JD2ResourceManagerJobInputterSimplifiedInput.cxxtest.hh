// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/resource_manager/ResourceOptionsFactory.cxxtest.hh
/// @brief test suite for basic::resource_manager::ResourceOptionsFactory
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <basic/resource_manager/LazyResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/Job.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <algorithm>
using namespace protocols::jd2;

class JD2ResourceManagerJobInputterSimplifiedInputTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	std::string input_v1(){
		return
			"Full:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <ResourceLocators>\n"
			"    <FileSystemResourceLocator tag=filesys/>\n"
			"  </ResourceLocators>\n"
			"  <ResourceOptions>\n"
			"    <SymmDataOptions tag=symm_data_options/>\n"
			"  </ResourceOptions>\n"
			"  <Resources>\n"
			"    <SymmData tag=1xu1FH_D.symm locator=filesys locatorID=1xu1FH_D.symm options=symm_data_options/>\n"
			"  </Resources>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata resource_tag=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}
	std::string input_v2(){
		return
			"wihtout options:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <ResourceLocators>\n"
			"    <FileSystemResourceLocator tag=filesys/>\n"
			"  </ResourceLocators>\n"
			"  <Resources>\n"
			"    <SymmData tag=1xu1FH_D.symm locator=filesys locatorID=1xu1FH_D.symm/>\n"
			"  </Resources>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata resource_tag=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}
	std::string input_v3(){
		return
			"Without ResourceLocator:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <Resources>\n"
			"    <SymmData tag=1xu1FH_D.symm locatorID=1xu1FH_D.symm/>\n"
			"  </Resources>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata resource_tag=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}
	std::string input_v4(){
		return
			"'file' rather than locatorID:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <Resources>\n"
			"    <SymmData tag=1xu1FH_D.symm file=1xu1FH_D.symm/>\n"
			"  </Resources>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata resource_tag=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}
	std::string input_v5(){
		return
			"getting tag from locatorID:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <Resources>\n"
			"    <SymmData file=1xu1FH_D.symm/>\n"
			"  </Resources>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata resource_tag=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}
	std::string input_v6(){
		return
			"moving resource definition into Data block:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <ResourceLocators>\n"
			"    <FileSystemResourceLocator tag=filesys/>\n"
			"  </ResourceLocators>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata resource_loader=SymmData locator=filesys locatorID=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}
	std::string input_v7(){
		return
			"moving resource definition into Data block shorter:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata resource_loader=SymmData file=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}
	std::string input_v8(){
		return
			"Figuring out resource loader from resource description:\n"
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=1xu1>\n"
			"      <Data desc=startstruct pdb=1xu1FH_D.pdb/>\n"
			"      <Data desc=symmdata file=1xu1FH_D.symm/>\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
	}

	void do_compare(
		std::string const & a,
		std::string const & b
	){
		basic::resource_manager::ResourceManager * resource_manager(
			basic::resource_manager::ResourceManager::get_instance());
		basic::resource_manager::LazyResourceManager * lazy_resource_manager(
			dynamic_cast< basic::resource_manager::LazyResourceManager * > (
			resource_manager ));

		// Make absolutely sure we're starting with a clean resource manager
		lazy_resource_manager->clear();

		std::istringstream a_jobstream(a);
		JD2ResourceManagerJobInputter a_inputter;
		Jobs a_jobvector;
		try {
			a_inputter.fill_jobs_from_stream( a_jobstream, a_jobvector );
			lazy_resource_manager->clear();
		} catch (utility::excn::Exception e ) {
			std::cerr << "The following resource manager XML definitions file produced an exception:" << std::endl;
			std::cerr << "XML a:" << std::endl << a << std::endl;
			std::cerr << "Raised exception:" << std::endl << e.msg() << std::endl;
			lazy_resource_manager->show(std::cerr);

			TS_ASSERT( false );
		}

		std::istringstream b_jobstream(b);
		JD2ResourceManagerJobInputter b_inputter;
		Jobs b_jobvector;
		try {
			b_inputter.fill_jobs_from_stream( b_jobstream, b_jobvector );
			lazy_resource_manager->clear();
		} catch (utility::excn::Exception e ) {
			std::cerr << "The following resource manager XML definitions file produced an exception:" << std::endl;
			std::cerr << "XML b:" << std::endl << b << std::endl;
			std::cerr << "Raised exception:" << std::endl << e.msg() << std::endl;
			lazy_resource_manager->show(std::cerr);

			TS_ASSERT( false );
			return;
		}

		for (
				Jobs::const_iterator
				ja = a_jobvector.begin(), jae = a_jobvector.end(), jb = b_jobvector.begin();
				ja != jae; ++ja, ++jb ) {
			if ( jb == b_jobvector.end() ) {
				std::cerr << "The following resource manager XML definitions produced different jobs:" << std::endl;
				std::cerr << "XML a:" << std::endl << a << std::endl;
				std::cerr << "XML b:" << std::endl << b << std::endl;

				std::cerr << "Job a: " << **ja << std::endl << std::endl;

				std::cerr << "Job b: NULL" << std::endl;
				TS_ASSERT( false );
				return;
			}
			if ( **ja != **jb ) {
				std::cerr << "The following resource manager XML definitions produced different jobs:" << std::endl;
				std::cerr << "XML a:" << std::endl << a << std::endl;
				std::cerr << "XML b:" << std::endl << b << std::endl;

				std::cerr << "Job a: " << **ja << std::endl << std::endl;
				std::cerr << "Job b: " << **jb << std::endl << std::endl;
				TS_ASSERT( false );
				return;
			}
		}
	}


	void test_compare_v1v2() {
		do_compare(input_v1(), input_v2());
	}

	void test_compare_v1v3() {
		do_compare(input_v1(), input_v3());
	}

	void test_compare_v1v4() {
		do_compare(input_v1(), input_v4());
	}

	void test_compare_v1v5() {
		do_compare(input_v1(), input_v5());
	}

	void test_compare_v1v6() {
		// TODO mjo: make this work
		//do_compare(input_v1(), input_v6());
	}

	void test_compare_v1v7() {
		// TODO mjo: make this work
		//do_compare(input_v1(), input_v7());
	}

	void test_compare_v1v8() {
		// TODO mjo: make this work
		//do_compare(input_v1(), input_v8());
	}
};
