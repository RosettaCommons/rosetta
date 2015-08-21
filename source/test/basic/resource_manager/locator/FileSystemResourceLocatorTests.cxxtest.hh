// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FileSystemResourceLocatorTests.cxxtest.hh
/// @brief  Test the FileSystemResourceLocator class
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>

#include <basic/Tracer.hh>

#include <basic/resource_manager/locator/FileSystemResourceLocator.hh>

// Utility headers
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

// C++ headers
#include <string>
#include <sstream>

static basic::Tracer tr( "basic.resource_manager.locator.FileSystemResourceLocator.cxxtest");

class FileSystemResourceLocatorTests : public CxxTest::TestSuite {

public:

	void test_FileSystemResourceLocator_locate_file() {
		using namespace basic::resource_manager;
		using namespace basic::resource_manager::locator;

		{
			utility::io::ozstream test_file("FileSystemResourceLocator.txt");
			test_file << "FileSystemResourceLocator test" << std::endl;
		}

		FileSystemResourceLocator file_locator;

		ResourceStreamOP rstream = file_locator.locate_resource_stream( "FileSystemResourceLocator.txt" );
		std::string test_contents;
		std::getline( rstream->stream(), test_contents );
		TS_ASSERT_EQUALS( test_contents, "FileSystemResourceLocator test" );

		TS_ASSERT( dynamic_cast< FileStream * > ( rstream.get() ) );
	}


};
