// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileLoaderCreator.hh
///
/// @brief      Loader class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)


#ifndef INCLUDED_protocols_membrane_io_LipoFileLoaderCreator_hh
#define INCLUDED_protocols_membrane_io_LipoFileLoaderCreator_hh

// Unit Headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>

namespace protocols {
namespace membrane {
namespace io {

/// @brief Lipo File Loader Creator
/// @details Register Lipo file loader with protocols init
class LipoFileLoaderCreator : public basic::resource_manager::ResourceLoaderCreator {

public:

    /// @brief Return lipo file loader class
	virtual
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const;

    /// @brief Return lipo file loader type
	virtual
	std::string
	loader_type() const;

}; // class LipoFileLoaderCreator

} // io
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_io_LipoFileLoaderCreator_hh

