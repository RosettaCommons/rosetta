// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <utility/SingletonBase.hh>

// C++ headers
#include <string>

namespace protocols {
namespace noesy_assign {

/// @brief WARNING WARNING WARNING THIS SINGLETON CLASS HOLDS NON-CONST
/// JOB-SPECIFIC DATA AND MAKES EVERYTHING THAT RELIES ON IT THREAD-UNSAFE.
/// THIS IS NOT HOW SINGLETONS SHOULD BE USED.
class CovalentCompliance : public utility::SingletonBase< CovalentCompliance >
{
public:
	friend class utility::SingletonBase< CovalentCompliance >;

private:
	/// @brief Private constructor for singleton class
	CovalentCompliance();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static CovalentCompliance * create_singleton_instance();

public:

	/// @brief This is clearly thread-unsafe.
	void load_dist_table( std::string const& file );
	bool is_compliant( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2, core::Real cutoff = 5.0 ) const;
	core::Real distance( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2 ) const;

private:
	FragsToAtomDistOP covalent_distances_;
};

}
}
#endif
