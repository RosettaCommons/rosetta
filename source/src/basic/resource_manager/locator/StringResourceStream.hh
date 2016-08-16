// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/StringResourceStream.hh
/// @brief  A thin wrapper around an std::stringstream for use with the ResourceManager
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_basic_resource_manager_locator_StringResourceStream_hh
#define INCLUDED_basic_resource_manager_locator_StringResourceStream_hh

//unit headers
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/locator/StringResourceStream.fwd.hh>

//utility headers

//C++ headers
#include <sstream>

namespace basic {
namespace resource_manager {
namespace locator {

/// @brief The %StringResourceStream is a wrapper class for a std::stringstream
/// that can be used to construct a resource.  Useful when reading the resource
/// stream in to memory (e.g. from a database or from multiple files) before
/// trying to construct the resource.
class StringResourceStream : public basic::resource_manager::ResourceStream
{
public:
	StringResourceStream();

	StringResourceStream(
		std::string const & contents
	);

	/// @brief Append a piece of data to the internal stringstream
	template< class InputSource >
	StringResourceStream(
		InputSource & in_stream
	) {
		in_stream >> stream_;
	}

private:
	/// @brief This is private. The StringResourceStream shouldn't be copied
	StringResourceStream(
		StringResourceStream const & src);

public:

	virtual
	~StringResourceStream();

	/// @brief Construct the stringstream from the given input string.
	virtual
	void
	fill(
		std::string const & contents);

	/// @brief Give non-const access to the internal stringstream so that it can be used to construct a resource.
	virtual
	std::istream &
	stream();

private: // members
	std::stringstream stream_;
};

} // namespace locator
} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_locator_StringResourceStream_hh
