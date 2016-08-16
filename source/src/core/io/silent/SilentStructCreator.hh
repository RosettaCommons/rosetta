// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentStructCreator.hh
/// @brief  Base class for SilentStructCreators for the SilentStruct load-time factory registration scheme
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_SilentStructCreator_hh
#define INCLUDED_core_io_silent_SilentStructCreator_hh

// Unit Headers
#include <core/io/silent/SilentStructCreator.fwd.hh>

// Package Headers
#include <core/io/silent/SilentStruct.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace core {
namespace io {
namespace silent {

/// @brief Abstract base class for a Mover factory; the Creator class is responsible for
/// creating a particular mover class.
class SilentStructCreator : public utility::pointer::ReferenceCount
{
public:
	SilentStructCreator();
	virtual ~SilentStructCreator();

	virtual SilentStructOP create_silent_struct() const = 0;
	virtual std::string keyname() const = 0;
};

} //namespace silent
} //namespace io
} //namespace core

#endif
