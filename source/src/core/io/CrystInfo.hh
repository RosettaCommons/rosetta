// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/CrystInfo.hh
/// @brief
/// @author  Frank DiMaio

#ifndef INCLUDED_core_io_CrystInfo_hh
#define INCLUDED_core_io_CrystInfo_hh

#include <core/types.hh>

#include <iostream>
#include <string>
#include <vector>

namespace core {
namespace io {

class CrystInfo
{
public:
	/// @brief default constructor to initialize all values
	CrystInfo() :
		A_( 0.0 ),
		B_( 0.0 ),
		C_( 0.0 ),
		alpha_( 90.0 ),
		beta_( 90.0 ),
		gamma_( 90.0 ),
		spacegroup_( "P 1" )
	{ }

	/// @brief Debug printing, serialazing to Tracer like object.
	friend std::ostream& operator <<(std::ostream &os, CrystInfo const & ci) {
		os << "<CrystInfo>{" << ci.A() << "," << ci.B() << "," << ci.C() << ","
			<< ci.alpha() << "," << ci.beta() << "," << ci.gamma() << " : " << ci.spacegroup() << "}";
		return os;
	}

	Real A() const { return A_; }
	void A(Real Ain) { A_=Ain; }
	Real B() const { return B_; }
	void B(Real Bin) { B_=Bin; }
	Real C() const { return C_; }
	void C(Real Cin) { C_=Cin; }

	Real alpha() const { return alpha_; }
	void alpha(Real alphain) { alpha_=alphain; }
	Real beta() const { return beta_; }
	void beta(Real betain) { beta_=betain; }
	Real gamma() const { return gamma_; }
	void gamma(Real gammain) { gamma_=gammain; }

	std::string spacegroup() const { return spacegroup_; }
	void spacegroup(std::string spacegroupin) { spacegroup_=spacegroupin; }


private:
	/// For now, all member names have the same names as fields in PDB standard.
	Real A_,B_,C_,alpha_,beta_,gamma_;
	std::string spacegroup_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace io
} // namespace core

#endif
