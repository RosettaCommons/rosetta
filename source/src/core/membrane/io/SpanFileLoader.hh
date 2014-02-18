// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileLoader.hh
///
/// @brief      Loader for SpanningTopology object - generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileLoader_hh
#define INCLUDED_core_membrane_io_SpanFileLoader_hh

// Unit Headers
#include <core/membrane/io/SpanFileLoader.fwd.hh>

// Project Headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <istream>

namespace core {
namespace membrane {
namespace io {
    
/// @brief Span File Loader
/// @details Load Spanning Topology Object
class SpanFileLoader : public basic::resource_manager::ResourceLoader {

public:

	/// @brief Constructor
	SpanFileLoader();

	/// @brief Destructor
	virtual ~SpanFileLoader();

	/// @brief Loads a SpanningTopology Object from spanfile
	virtual
	utility::pointer::ReferenceCountOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
		) const;

	/// @brief Default Options
	virtual
	basic::resource_manager::ResourceOptionsOP
	default_options() const;

}; // class SpanFileLoader

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_SpanFileLoader_hh

