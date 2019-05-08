// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/SharedSilentData.hh
///
/// @brief silent input file reader for mini
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_SharedSilentData_hh
#define INCLUDED_core_io_silent_SharedSilentData_hh

// mini headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
}; // class SharedSilentData

class SimpleSequenceData : public SharedSilentData {

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
}; // class SimpleSequenceData


} // namespace silent
} // namespace io
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_io_silent_SharedSilentData )
#endif // SERIALIZATION

#endif
