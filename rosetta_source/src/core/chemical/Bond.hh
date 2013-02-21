// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin Bond
///
/// @brief
/// A class for defining chemical Bonds, with properties specific to a ResidueType, not conformation info
/// specific to a Residue. Conformation info goes in conformation::Bond. BondTypes are not ResidueType specific.
///
///
///
///
/// @author
/// Gordon Lemmon
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_Bond_hh
#define INCLUDED_core_chemical_Bond_hh


// Unit headers
#include <core/chemical/Bond.fwd.hh>
#include <core/types.hh>

namespace core {
namespace chemical {

/// @brief basic chemical Bond
///
/// @details name, element, certain properties and parameters from .params file
///
class Bond {

public:
	Bond(): distance_(0), bond_name_(BondName()){}
	Bond(Real d, BondName name): distance_(d), bond_name_(name){}

// Setters
	void distance(Real distance){
		distance_ = distance;
	}
	void bond_name(BondName bond_name){
		bond_name_ = bond_name;
	}

// Setters
	Real distance(){
		return distance_;
	}
	BondName bond_name(){
		return bond_name_;
	}

	void print( std::ostream & out ) const;

	friend
	std::ostream & operator<<(std::ostream & out, Bond const & bond );

private:

	Real distance_;
	BondName bond_name_; // this is an enum defined in ResidueType.fwd.hh
};


} // chemical
} // core



#endif // INCLUDED_core_chemical_Bond_HH
