// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/database_io.cc
/// @brief   Database input/output function definitions for carbohydrate-specific data.
/// @author  Labonte

// Unit header
#include <core/chemical/carbohydrates/database_io.hh>

// Project headers
#include <core/types.hh>

// Utility headers
//#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>

// Basic headers
#include <basic/Tracer.hh>


// Construct tracer.
static basic::Tracer TR("core.chemical.carbohydrates.database_io");


namespace core {
namespace chemical {
namespace carbohydrates {

using namespace std;
using namespace utility;
using namespace utility::file;
using namespace utility::io;
using namespace core;

// Return a list of strings, which are saccharide-specific properties and modifications, read from a database file.
utility::vector1<std::string>
read_properties_from_database_file(std::string const & filename)
{
	// Check if file exists.
	if(!file_exists(filename)) {
		utility_exit_with_message("Cannot find database file: '" + filename + "'");
	}

	// Open file.
	izstream data(filename.c_str());
	if (!data.good()) {
		utility_exit_with_message("Unable to open database file: '" + filename + "'");
	}

	string line;
	vector1<string> properties;

	// Read lines.
	while (getline(data, line)) {
		trim(line, " \t\n");  // Remove leading and trailing whitespace.
		if ((line.size() < 1) || (line[0] == '#')) continue;  // Skip comments and blank lines.

		properties.push_back(line);
	}

	// Close file.
	data.close();

	TR.Debug << "Read " << properties.size() << " properties from the carbohydrate database." << endl;

	return properties;
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
