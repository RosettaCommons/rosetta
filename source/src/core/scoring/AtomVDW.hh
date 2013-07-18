// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/VDW_Energy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_AtomVDW_hh
#define INCLUDED_core_scoring_AtomVDW_hh

// Unit Headers
#include <core/scoring/AtomVDW.fwd.hh>

// Package headers
#include <core/chemical/AtomTypeSet.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

namespace core {
namespace scoring {


class AtomVDW : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor, reads data file
	AtomVDW( std::string const & atom_type_set_name );


	///
	utility::vector1< Real > const &
	operator()( Size const atom_type_index ) const
	{
		return atom_vdw_[ atom_type_index ];
	}

	///
	Real
	approximate_vdw_radius( Size const atom_type_index ) const
	{
		return approximate_vdw_radii_[ atom_type_index ];
	}

	std::string atom_type_set_name() const;

private:

	void
	setup_approximate_vdw_radii(
		utility::vector1< int > const & atom_type_index,
		chemical::AtomTypeSet const & atom_type_set
	);

private:
	std::string atom_type_set_name_;
	utility::vector1< utility::vector1< Real > > atom_vdw_;

	/// @brief  Approximation to per-atom radii, derived from atom_vdw_ data
	utility::vector1< Real > approximate_vdw_radii_;

};

}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
