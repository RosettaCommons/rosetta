// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileOptions.hh
///
/// @brief      Options class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_io_LipoFileOptions_hh
#define INCLUDED_protocols_membrane_io_LipoFileOptions_hh

// Unit Headers
#include <protocols/membrane/io/LipoFileOptions.fwd.hh>

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

/// @brief Lipo File Options
/// @details Options for reading in lipid acc data
class LipoFileOptions : public basic::resource_manager::ResourceOptions {

public:

    /// @brief Constructor
	LipoFileOptions();
    
    /// @brief Destructor
	virtual ~LipoFileOptions();

    /// @brief Parse options from .xml definition file
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
		);

    /// @brief Return options class type
	virtual
	std::string
	type() const;

}; // class LipoFileOptions

} // io
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_io_LipoFileOptions_hh
