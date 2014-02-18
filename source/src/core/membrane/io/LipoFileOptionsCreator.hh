// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileOptionsCreator.hh
///
/// @brief      Options creator class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileOptionsCreator_hh
#define INCLUDED_core_membrane_io_LipoFileOptionsCreator_hh

//package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Options Creator for Lipo Files
/// @details Register Lipo File Options with core init
class LipoFileOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{

public:

    /// @brief Return lipo file options type
	virtual std::string options_type() const;
    
    /// @brief Return lipo file options class
	virtual basic::resource_manager::ResourceOptionsOP create_options() const;

}; // class LipoFileOptionsCreator

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileOptionsCreator_hh

