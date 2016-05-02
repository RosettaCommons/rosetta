// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/resource_manager/JD2ResourceManagerJobInputter.cxxtest.hh
/// @brief test suite for protocols::jd2::JD2ResourceManagerJobInputter
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/Job.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

using namespace protocols::jd2;


class JD2ResourceManagerJobInputterTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_JD2ResourceManagerJobInputter_read_job_w_nstruct_10() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb />\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TS_ASSERT( jobvector.size() == 10 );
	}

	void test_JD2ResourceManagerJobInputter_read_job_w_nstruct_10_option_without_namespace() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb />\n"
			"      <Option native=1ten.pdb />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TS_ASSERT( jobvector.size() == 10 );
	}

	void test_JD2ResourceManagerJobInputter_missing_desc_in_job() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data pdb=1ten.pdb />\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Failed to find a data description (desc) amongst the options pairs listed reading a 'Data' tag in a Job tag.\n A desc option must always be given.\n"
				"Error:  'pdb' tag given for a non-'startstruct' option in the 'Data' tag of a Job tag.\nProblem encountered for job named 'firstjob'.\n"
				"Options given:\n"
				"\t(pdb, 1ten.pdb)\n"
				"Thrown from protocols::jd2::JD2ResourceManagerJobInputter::parse_job_tag\n";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}

	}

	void test_JD2ResourceManagerJobInputter_missing_resource_and_pdb_for_startstruct_data_in_job() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct />\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Failed to find a resource name or a pdb name reading a 'Data' tag in a Job tag.  Either a 'resource_tag' or a 'pdb' option must be provided.\n"
				"Problem encountered for job named 'firstjob'.\n"
				"Options given:\n"
				"\t(desc, startstruct)\n"
				"Thrown from protocols::jd2::JD2ResourceManagerJobInputter::parse_job_tag\n";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManagerJobInputter_found_both_resource_and_pdb_for_startstruct_data_in_job() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb resource_tag=1ten_native/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Error: Both a 'resource_tag' and a 'pdb' tag were found for a 'Data' tag in the Job tag.\nProblem encountered for job named 'firstjob'.\n"
				"Options given:\n"
				"\t(desc, startstruct)\n"
				"\t(pdb, 1ten.pdb)\n"
				"\t(resource_tag, 1ten_native)\n"
				"Thrown from protocols::jd2::JD2ResourceManagerJobInputter::parse_job_tag\n";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error );
		}
	}

	/// Only the starstruct is allowed to give a "pdb" option
	void test_JD2ResourceManagerJobInputter_pdb_without_startstruct_in_job() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=whatever pdb=1ten.pdb />\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Error:  'pdb' tag given for a non-'startstruct' option in the 'Data' tag of a Job tag.\nProblem encountered for job named 'firstjob'.\n"
				"Options given:\n"
				"\t(desc, whatever)\n"
				"\t(pdb, 1ten.pdb)\n"
				"Thrown from protocols::jd2::JD2ResourceManagerJobInputter::parse_job_tag\n";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error );
		}

	}

	void test_JD2ResourceManagerJobInputter_job_without_startstruct() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=first nstruct=10 >\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error = "Error: Job 'first' given without a 'startstruct'";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error );
		}

	}

	void test_JD2ResourceManagerJobInputter_integer_option_read() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option packing:ex1:level=4 />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cout << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TS_ASSERT( jobvector.size() == 10 );
	}

	void test_JD2ResourceManagerJobInputter_integer_option_read_error() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option packing:ex1:level=four />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Error converting value 'four' given for option 'packing:ex1:level' to an integer from within JD2ResourceManagerJobInputter::parse_job_tag\n";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error );
		}

	}
	void test_JD2ResourceManagerJobInputter_real_option_read() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option out:scorecut=0.3 />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cout << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TS_ASSERT( jobvector.size() == 10 );

	}
	void test_JD2ResourceManagerJobInputter_real_option_read_error() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option out:scorecut=0fge.3 />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Error converting value '0fge.3' given for option 'out:scorecut' to a floating point number from within JD2ResourceManagerJobInputter::parse_job_tag\n";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error );
		}

	}
	void test_JD2ResourceManagerJobInputter_boolean_option_read() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option packing:ex1=1 />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			TS_ASSERT( false );
			std::cout << e.msg() << std::endl;
		}
		TS_ASSERT( jobvector.size() == 10 );

	}
	void test_JD2ResourceManagerJobInputter_boolean_option_read_error() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option packing:ex1=true />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Error converting value 'true' given for option 'packing:ex1' to a boolean from within JD2ResourceManagerJobInputter::parse_job_tag\n"
				" Boolean options must be given either a '1' or a '0'";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}

			TS_ASSERT( e.msg() == expected_error );
		}

	}

	void test_JD2ResourceManagerJobInputter_nonexistant_option_read_error() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option packing::ex1=1 />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
			TS_ASSERT( false ); // should throw an exception
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Error: Option key 'packing::ex1' not found. Please remember to use only one colon when giving options.\n";
			if ( e.msg() != expected_error ) {
				std::cout << "expected error: '" << expected_error << "'" << std::endl;
				std::cout << "actual error:   '" << e.msg() << "'" << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}

	}

	void test_JD2ResourceManagerJobInputter_nonexistant_option_read_error2() {
		JD2ResourceManager * jd2rm(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		jd2rm->clear();

		std::string xmlfile =
			"<JD2ResourceManagerJobInputter>\n"
			"  <Jobs>\n"
			"    <Job name=firstjob nstruct=10 >\n"
			"      <Data desc=startstruct pdb=1ten.pdb/>\n"
			"      <Option in:file:native=1ten.pdb />\n"
			"      <Option ex1=1 />\n"
			"    </Job>\n"
			"  </Jobs>\n"
			"</JD2ResourceManagerJobInputter>\n";
		std::istringstream jobstream( xmlfile );
		JD2ResourceManagerJobInputter inputter;
		Jobs jobvector;
		try {
			inputter.fill_jobs_from_stream( jobstream, jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			TS_ASSERT( false ); // shouldn't throw an exception
		}

	}

};
