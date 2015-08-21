// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for holding bond information.
///
/// @author
/// Gordon Lemmon, Rocco Moretti (rmorettiase@gmail.com)
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_Bond_hh
#define INCLUDED_core_chemical_Bond_hh


// Unit headers
#include <core/chemical/Bond.fwd.hh>
#include <core/types.hh>

#include <string>

namespace core {
namespace chemical {

BondName convert_to_BondName( std::string const & id );

/// @brief basic chemical Bond
///
/// @details name, element, certain properties and parameters from .params file
///
class Bond {

public:

	Bond();
	Bond(Real d, BondName name, bool cut_bond = false);
	Bond(Real d, BondOrder order, BondConjugability conj,  BondAromaticity aroma, BondRingness ring = UnknownRingness, BondIsometry isom = UnknownIsometry, bool cut_bond = false);

	// Setters
	void distance(Real distance){
		distance_ = distance;
	}
	void bond_name(BondName bond_name);

	void cut_bond(bool cut_bond){
		cut_bond_ = cut_bond;
	}
	void order(BondOrder order);

	void conjugability(BondConjugability conjug){
		conjug_ = conjug;
	}
	void aromaticity(BondAromaticity aroma);

	void ringness(BondRingness ring){
		ring_ = ring;
	}
	void isometry(BondIsometry isom){
		isometry_ = isom;
	}


	// Getters
	Real distance() const {
		return distance_;
	}
	BondName bond_name() const {
		return bond_name_;
	}
	bool cut_bond() const {
		return cut_bond_;
	}
	BondOrder order() const {
		return order_ ;
	}
	BondConjugability conjugability() const {
		return conjug_ ;
	}
	BondAromaticity aromaticity() const {
		return aroma_ ;
	}
	BondRingness ringness() const {
		return ring_ ;
	}
	BondIsometry isometry() const {
		return isometry_ ;
	}

	// General functions

	void print( std::ostream & out ) const;

	friend
	std::ostream & operator<<(std::ostream & out, Bond const & bond );

	// Derived Data Access.

	/// @brief Reset the internal data such that it matches the appropriate value for the SDF datatype.
	///
	/// Substitution (taken from) for BCL's FindBondTypeFromSDFInfo( const std::size_t &SDF_ID)
	void SetSDFType( const core::Size SDF_ID);

	core::Size GetNumberOfElectrons() const;
	core::Size GetMinimumElectrons() const;
	core::Size GetMaximumElectrons() const;

	bool IsBondOrderKnown() const { return order_ != UnknownBondOrder; }

	bool IsBondInRing() const { return ring_ == BondInRing; }

	core::Size GetSDFileID() const;
	core::Size GetSDAltFileID() const;

	/// @brief Return true if this bond represents a non-physical bond
	bool is_fake() const {
		return order_ == OrbitalBondOrder || order_ == PseudoBondOrder;
		// An UnknownBondOrder bond is treated as a real bond (of unknown order).
	}

private:

	Real distance_;
	/// @brief Is the bond allowed to be in the atom tree?
	bool cut_bond_;
	// The following enums are defined in Bond.fwd.hh
	BondName bond_name_;
	/// @brief The bond order (single double triple ..)
	BondOrder order_;
	/// @brief Can the bond participate in a conjugated system
	BondConjugability conjug_;
	/// @brief Is the bond in an aromatic ring? (Distinct from being conjugatable and in a ring.)
	BondAromaticity aroma_;
	/// @brief Is the bond in a ring?
	BondRingness ring_;
	/// @brief For double bonds, what's the E/Z isometry?
	BondIsometry isometry_;
};


} // chemical
} // core


#endif // INCLUDED_core_chemical_Bond_HH
