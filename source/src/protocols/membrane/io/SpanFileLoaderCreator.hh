// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileLoaderCreator.hh
///
/// @brief      Loader for SpanningTopology object - generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_io_SpanFileLoaderCreator_hh
#define INCLUDED_protocols_membrane_io_SpanFileLoaderCreator_hh

// Unit Headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>

namespace protocols {
namespace membrane {
namespace io {

/// @brief Span File Loader Creator Class
/// @details Register span file loader with protocols init
class SpanFileLoaderCreator : public basic::resource_manager::ResourceLoaderCreator {

public:

    /// @brief Return span file resource loader class
	virtual
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const;

    /// @brief Return loader type
	virtual
	std::string
	loader_type() const;

}; // class SpanFileLoader

} // io
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_io_SpanFileLoaderCreator_hh

