// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/locator/StringResourceStream.hh
/// @brief  A thin wrapper around an std::stringstream for use with the ResourceManager
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_basic_resource_manager_locator_StringResourceStream_hh
#define INCLUDED_basic_resource_manager_locator_StringResourceStream_hh

//unit headers
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/locator/StringResourceStream.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>

//C++ headers
#include <sstream>

namespace basic {
namespace resource_manager {
namespace locator {

class StringResourceStream : public basic::resource_manager::ResourceStream
{
public:
	StringResourceStream();

	StringResourceStream(
		std::string const & contents
	);

	template< class InputSource >
	StringResourceStream(
		InputSource & in_stream
	) {
		in_stream >> stream_;
	}

private:
	///@brief This is private. The StringResourceStream shouldn't be copied
	StringResourceStream(
		StringResourceStream const & src);

public:

	virtual
	~StringResourceStream();

	virtual
	void
	fill(
		std::string const & contents);

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
