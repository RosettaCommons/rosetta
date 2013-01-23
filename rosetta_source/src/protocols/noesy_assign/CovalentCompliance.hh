// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/noesy_assign/CovalentCompliance.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_CovalentCompliance_HH
#define INCLUDED_protocols_noesy_assign_CovalentCompliance_HH


// Unit Headers
//#include <protocols/noesy_assign/CovalentCompliance.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/FragsToAtomDist.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>

// Utility headers

// C++ headers
#include <string>

namespace protocols {
namespace noesy_assign {

class CovalentCompliance {
private:
  //Singleton Class
  CovalentCompliance();

public:
  static CovalentCompliance const* get_instance();
  static CovalentCompliance * get_nonconst_instance();
  void load_dist_table( std::string const& file );
  bool is_compliant( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2, core::Real cutoff = 5.0 ) const;
  core::Real distance( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2 ) const;

private:
  static CovalentCompliance* instance_;
  FragsToAtomDistOP covalent_distances_;
};

}
}
#endif
