// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileOptionsCreator.hh
///
/// @brief      Options for generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_io_SpanFileOptionsCreator_hh
#define INCLUDED_protocols_membrane_io_SpanFileOptionsCreator_hh

//package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace membrane {
namespace io {

/// @brief Options for generating per-chain membrane spanning data
/// @details Generate spanning topology object options set
class SpanFileOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{

public:

    /// @brief Return options class type
	virtual std::string options_type() const;
    
    /// @brief Return options class
	virtual basic::resource_manager::ResourceOptionsOP create_options() const;

}; // class SpanFileOptionsCreator

} // io
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_io_SpanFileOptionsCreator_hh

