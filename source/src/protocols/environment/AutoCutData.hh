// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/AutoCutData.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_AutoCutData_hh
#define INCLUDED_protocols_environment_AutoCutData_hh

// Unit Headers
#include <protocols/environment/AutoCutData.fwd.hh>
#include <protocols/environment/AutoCutDataCreator.hh>

// Package headers

// Project headers
#include <basic/datacache/WriteableCacheableData.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <set>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class AutoCutData : public basic::datacache::WriteableCacheableData {
  typedef basic::datacache::WriteableCacheableData Parent;
  typedef basic::datacache::WriteableCacheableDataOP ParentOP;

public:
  AutoCutData( std::istream &in );

  AutoCutData( core::Size const& hash,
               std::set< Size > const& cuts );

  void write( std::ostream &out ) const;

  basic::datacache::CacheableDataOP clone() const;

  std::string datatype() const;

  utility::vector1< Size > const& cuts() const{
    return auto_cuts_;
  }

  core::Size const& hash() const {
    return hash_;
  }

private:
  utility::vector1< core::Size > auto_cuts_;
  core::Size hash_;
};

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_AutoCutData_hh
