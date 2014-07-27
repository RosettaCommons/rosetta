// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/LocalPosition.hh
/// @brief A class that is used to express a mover-specific DoF-unlock on a ProtectedPose. Its destruction expresses a re-locking.
///
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_LocalPosition_hh
#define INCLUDED_protocols_environment_LocalPosition_hh

// Unit Headers
#include <core/environment/LocalPosition.fwd.hh>

// Package headers

// Project headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

// C++ Headers
#include <string>

// ObjexxFCL Headers

namespace core {
namespace environment {

class LocalPosition : public utility::pointer::ReferenceCount {

public:
  LocalPosition();

  LocalPosition( std::string const& comma_deliniated );

  LocalPosition( std::string const&, core::Size const& );

  std::string const& label() const;

  core::Size const& position() const;

  bool operator<( LocalPosition const& ) const;

  bool operator==( LocalPosition const& ) const;

  bool operator!=( LocalPosition const& ) const;

//  bool operator>( LocalPosition const& ) const;
//
//  bool operator>=( LocalPosition const& ) const;
//
//  bool operator<=( LocalPosition const& ) const;

protected:
  void label( std::string const& );

  void position( core::Size );

private:
  std::string label_;

  core::Size position_;

}; // end LocalPosition base class

static const LocalPosition NO_POSITION = LocalPosition( "", 0 );

extern std::ostream& operator<<( std::ostream&, LocalPosition const& );

} // environment
} // core

#endif //INCLUDED_protocols_environment_LocalPosition_hh
