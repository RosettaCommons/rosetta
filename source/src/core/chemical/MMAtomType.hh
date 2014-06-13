// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/MMAtomType.hh
/// @brief  Molecular mechanics atom type class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_chemical_MMAtomType_hh
#define INCLUDED_core_chemical_MMAtomType_hh

// Unit headers
#include <core/chemical/MMAtomType.fwd.hh>

// Project headers
#include <core/types.hh>

// C++ headers
#include <string>

namespace core {
namespace chemical {

/// @brief Basic MM atom type
///
/// @details Simple class for holding the name and the LJ properties of a Charmm
/// molecular mechanics atom type. Borrows heavily and functions similarly
/// to the rosetta atom type class, AtomType
///
class MMAtomType
{

public:

	///  @brief Construct a new MMAtomType with its name
	MMAtomType( std::string const & name_in):
		name_( name_in ),
		lj_radius_( 0.0 ),
		lj_wdepth_( 0.0 ),
		lj_three_bond_radius_( 0.0 ),
		lj_three_bond_wdepth_( 0.0 )
	{}

	/// @brief Return the name of the MMAtomType
	std::string const& name() const { return name_; }

	/// @brief Return the LJ radius of the atom type
	Real lj_radius() const { return lj_radius_; }

	/// @brief Return the LJ well depth of the atom type
	Real lj_wdepth() const { return lj_wdepth_; }

	/// @brief Return the LJ radius for use when atoms types are seperated by 3 bonds
	Real lj_three_bond_radius() const { return lj_three_bond_radius_; }

	/// @brief Return the LJ well depth for use when atoms types are seperated by 3 bonds
	Real lj_three_bond_wdepth() const { return lj_three_bond_wdepth_; }

	/// @brief set LJ and LK solvation parameter for this atom type
	void set_parameter(	std::string const & param, Real const setting	);

private:

	/// @brief name of the mm atom type
	std::string const name_;

	/// @brief Lennard-Jones parameters
	Real lj_radius_;
	Real lj_wdepth_;
	Real lj_three_bond_radius_;
	Real lj_three_bond_wdepth_;

};


} // chemical
} // core



#endif // INCLUDED_core_chemical_MMAtomType_HH
