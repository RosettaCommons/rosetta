// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/deallocation/ResourceDeallocationMessage.hh
/// @brief  The definition for class protocols::jd3::ResourceDeallocationMessage
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_deallocation_ResourceDeallocationMessage_HH
#define INCLUDED_protocols_jd3_deallocation_ResourceDeallocationMessage_HH

// Unit headers
#include <protocols/jd3/deallocation/ResourceDeallocationMessage.fwd.hh>

// Package headers
#include <protocols/jd3/deallocation/DeallocationMessage.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace deallocation {

/// @brief %ResourceDeallocationMessage class is used so that a JobQueen that knows a resource will
/// no longer be needed can communicate this knowledge to other JobQueens running remotely so that
/// they can deallocate that resource.
class ResourceDeallocationMessage : public DeallocationMessage
{
public:

	ResourceDeallocationMessage();
	ResourceDeallocationMessage( std::string const & resource_name );
	virtual ~ResourceDeallocationMessage();

	std::string const & resource_name() const;
	void resource_name( std::string const & setting );

private:
	std::string resource_name_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // ResourceDeallocationMessage

} // namespace deallocation
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_deallocation_ResourceDeallocationMessage )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_deallocation_ResourceDeallocationMessage_HH
