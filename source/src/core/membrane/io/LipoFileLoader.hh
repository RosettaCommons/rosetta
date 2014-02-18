// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileLoader.hh
///
/// @brief      Loader class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)


#ifndef INCLUDED_core_membrane_io_LipoFileLoader_hh
#define INCLUDED_core_membrane_io_LipoFileLoader_hh

// Unit Header
#include <core/membrane/io/LipoFileLoader.fwd.hh>

// Project Headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <istream>

namespace core {
namespace membrane {
namespace io {

/// @brief Load Lipid Acc data object
/// @details Load lipid acc data based on spanning topology from .lips4
class LipoFileLoader : public basic::resource_manager::ResourceLoader {

public:

	/// @brief Constructor
	LipoFileLoader();

	/// @brief Destructor
	virtual ~LipoFileLoader();

	/// @brief Returns an OP to lipo file loader options
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

}; // class LipoFileLoader

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileLoader_hh

