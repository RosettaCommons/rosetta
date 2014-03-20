// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileIO.hh
///
/// @brief      Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileIO_hh
#define INCLUDED_core_membrane_io_LipoFileIO_hh

// Unit Headers
#include <core/membrane/io/LipoFileIO.fwd.hh>

// Project Headers
#include <core/membrane/properties/LipidAccInfo.hh>

// Package Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/tag/Tag.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace core::membrane::properties;

namespace core {
namespace membrane {
namespace io {

/// @brief Lipid Acc Data IO
/// @details Read lipid accessibility data from .lips4 file
class LipoFileIO : public utility::pointer::ReferenceCount {

public:

    /// @brief Constructor
	LipoFileIO();

    /// @brief Destructor
	~LipoFileIO();

	/// @brief Copy Constructor
	/// @param Lipid Builder Object
	// Undefined, commenting out to fix PyRosetta build  LipoFileIO( LipoFileIO const & src );

	/// @brief Main lipid object construction function
	/// @param [ lipofile ]
	LipidAccInfoOP get_lips_exp_from_lipofile( std::string lipofile );

private: // functions

	/// @brief Read and store data from lips file
	/// @param Lipid Info object, lipid info file, and pose
	void read_lips_data(
		LipidAccInfoOP lipid_exp,
		std::string lipsfile
		);

}; // class LipoFileIO

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileIO_hh
