// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileIO.hh
///
/// @brief      Generates per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileIO_hh
#define INCLUDED_core_membrane_io_SpanFileIO_hh

// Unit Headers
#include <core/membrane/io/SpanFileIO.fwd.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/types.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

using namespace core::conformation::membrane;

namespace core {
namespace membrane {
namespace io {

/// @brief   Span File Reader
/// @details Reads in topology data from an octopus span file
///          for helical transmembrane proteins
class SpanFileIO : public utility::pointer::ReferenceCount {

public:

	/// @brief Constructor
	SpanFileIO();

    /// @brief Destructor
	~SpanFileIO();

	/// @brief Copy constructor
	// Undefined, commenting out to fix PyRosetta build  SpanFileIO( SpanFileIO const & src );

	/// @brief Get Topology from Spanfile
	/// @param spanfile
	SpanningTopologyOP get_topology_from_spanfile( std::string spanfile );

private: // functions

	/// @brief Read in information from given file
	/// @param Topology Info Object and spanfile
	void read_spanfile( SpanningTopologyOP topology, std::string spanfile );

	/// @brief Set up span and full span info
	/// @param Topology Info Object
	void setup_span_info( SpanningTopologyOP topology );

	/// @brief Set values in tmregion (true/false for is_tmregion)
	/// @param Topology Info Object
	void setup_tmregion( SpanningTopologyOP topology );

	/// @brief Set relative tmh information
	/// @param Topology Info Object
	void setup_relative_tmh( SpanningTopologyOP topology );

}; // class SpanFileIO

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_SpanFileIO_hh
