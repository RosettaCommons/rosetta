// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/Raw.fwd.hh
///
/// @brief Forward declarations for raw data classes
/// @author James Thompson, Monica Berrondo

#ifndef INCLUDED_core_io_raw_data_Raw_fwd_hh
#define INCLUDED_core_io_raw_data_Raw_fwd_hh

// mini headers
// AUTO-REMOVED #include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>

// ObjexxFCL Headers

// C++ Headers
#include <map>

namespace core {
namespace io {
namespace raw_data {

	// abstract base classes
	class RawStruct;
	class RawFileData;

	// derived classes
	class DecoyStruct;
	class DecoyFileData;
	class ScoreStruct;
	class ScoreFileData;

	// owning pointers
	typedef utility::pointer::shared_ptr< RawStruct > RawStructOP;
	typedef utility::pointer::shared_ptr< DecoyStruct > DecoyStructOP;
	typedef utility::pointer::shared_ptr< ScoreStruct > ScoreStructOP;

	// data types
	typedef std::map< std::string, RawStructOP > StructureMap;
} // namespace silent
} // namespace io
} // namespace core

#endif
