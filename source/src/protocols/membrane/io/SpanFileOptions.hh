// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileOptions.hh
///
/// @brief      Options for generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_io_SpanFileOptions_hh
#define INCLUDED_protocols_membrane_io_SpanFileOptions_hh

// Unit Headers
#include <protocols/membrane/io/SpanFileOptions.fwd.hh>

// Project Headers
#include <basic/resource_manager/ResourceOptions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <cstdlib>

namespace protocols {
namespace membrane {
namespace io {

/// @brief Options for generating per-chain membrane spanning data
/// @details Generate spanning topology object options set
class SpanFileOptions : public basic::resource_manager::ResourceOptions {

public:
    
    /// @brief Constructor
	SpanFileOptions();
    
    /// @brief Destructor
	virtual ~SpanFileOptions();

    /// @brief Parse options from .xml definition
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
		);

    /// @brief Return options class type
	virtual
	std::string
	type() const;

}; // class SpanFileOptions

} // io
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_io_SpanFileOptions_hh
