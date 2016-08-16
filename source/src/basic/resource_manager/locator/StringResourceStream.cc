// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/StringResourceStream.cc
/// @brief  A thin wrapper around an std::stringstream for use with the ResourceManager
/// @author Matthew O'Meara (mattjomeara@gmail.com)

//unit headers
#include <basic/resource_manager/locator/StringResourceStream.hh>

//C++ Headers
#include <sstream>

namespace basic {
namespace resource_manager {
namespace locator {

using std::string;
using std::stringstream;
using std::istream;

/// @detail This is private. The StringResourceStream shouldn't be copied
StringResourceStream::StringResourceStream(
	StringResourceStream const & src
) :
	ResourceStream( src ),
	stream_()
{}

/// @detail If you use this constructor be sure to set the string
///before accessing the stream
StringResourceStream::StringResourceStream() :
	ResourceStream(),
	stream_()
{}

StringResourceStream::StringResourceStream(
	string const & contents
) :
	stream_(contents)
{}

StringResourceStream::~StringResourceStream() {}

void
StringResourceStream::fill(
	string const & contents
) {
	stream_ << contents;
}

istream &
StringResourceStream::stream() {
	return stream_;
}


} // namespace locator
} // namespace resource_manager
} // namespace basic
