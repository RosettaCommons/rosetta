// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/environment/EnvCore.hh
/// @brief the low compile-level functions of Environment.
/// @details This class is here so that DofPassport can be in core and have a private constructor, while still satisfying the compile-order rules (Environment needs to be able to talk to movers, which are in protocols.1). It is AGAINST THE LAW to inherit directly from this class. Instead, inheret from Environment, if you really want to make a special kind of environment.
/// @author Justin Porter

#ifndef INCLUDED_core_environment_EnvCore_hh
#define INCLUDED_core_environment_EnvCore_hh

// Unit Headers
#include <core/environment/EnvCore.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// Package headers
#include <core/environment/DofPassport.fwd.hh>

// Project headers
#include <core/types.hh>

// C++ Headers
#include <ostream>

// ObjexxFCL Headers

namespace core {
namespace environment {

class EnvCore : public utility::pointer::ReferenceCount {
public:
	EnvCore( std::string const& env_name );

	virtual ~EnvCore();

	std::string const& name() const;

	EnvCoreCAP superenv() const;

	core::Size const& id() const;

protected:

	DofPassportOP issue_passport( std::string const& mover_name ) const;

	void set_superenv( EnvCoreCAP );

private:

	static core::Size generate_id();

	std::string const name_;

	core::Size id_;

	EnvCoreCAP superenv_;

	//TODO: make this thread safe
	static core::Size current_maximum_id_; // = 0

}; // end EnvCore base class

} // environment
} // core

#endif //INCLUDED_core_environment_EnvCore_HH
