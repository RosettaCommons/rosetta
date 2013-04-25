// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopsFileLoader.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_loops_LoopsFileLoader_HH
#define INCLUDED_protocols_loops_LoopsFileLoader_HH

//unit headers
#include <protocols/loops/LoopsFileLoader.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>


//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>

namespace protocols {
namespace loops {

class LoopsFileLoader : public basic::resource_manager::ResourceLoader
{
public:
	LoopsFileLoader();
	virtual ~LoopsFileLoader();

	/// @brief Returns a pointer to a WrappedPrimative< LoopsFileData > object
	/// which is constructed from the given input stream (istream) which in tern
	/// originates from a particular data source (given by the name input_tag)
	virtual
	utility::pointer::ReferenceCountOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
	) const;

	virtual
	basic::resource_manager::ResourceOptionsOP
	default_options() const;

};

} // namespace loops
} // namespace protocols



#endif //INCLUDED_protocols_loops_LoopsFileLoader_HH
