// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/SharedSilentData.hh
///
/// @brief silent input file reader for mini
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_SharedSilentData_hh
#define INCLUDED_core_io_silent_SharedSilentData_hh

// mini headers
// AUTO-REMOVED #include <core/types.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

namespace core {
namespace io {
namespace silent {

enum SharedSilentDataType {
	energynames = 1,
	simplesequencedata
};

class SharedSilentData : public utility::pointer::ReferenceCount {

public:
	SharedSilentData() {}

}; // class SharedSilentData

class SimpleSequenceData : public	SharedSilentData {

public:

	SimpleSequenceData() :
		sequence_( "" )
	{}

	void set_sequence( std::string sequence ) {
		sequence_ = sequence;
	}

	std::string sequence() {
		return sequence_;
	}

private:
	std::string sequence_;
}; // class SimpleSequenceData


} // namespace silent
} // namespace io
} // namespace core

#endif
