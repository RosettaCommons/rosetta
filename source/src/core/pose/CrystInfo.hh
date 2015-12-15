// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/CrystInfo.hh
/// @brief
/// @author  Frank DiMaio

#ifndef INCLUDED_core_pose_CrystInfo_hh
#define INCLUDED_core_pose_CrystInfo_hh

#include <core/types.hh>

#include <iostream>
#include <string>
#include <vector>

namespace core {
namespace pose {

class CrystInfo
{
public:
	/// @brief default constructor to initialize all values
	CrystInfo() {
		A_ = B_ = C_ = 0.0;
		alpha_ = beta_ = gamma_ = 90.0;
		spacegroup_ = "P 1";
	}

	/// @brief Debug printing, serialazing to Tracer like object.
	friend std::ostream& operator <<(std::ostream &os, CrystInfo const & ci) {
		os << "<CrystInfo>{" << ci.A() << "," << ci.B() << "," << ci.C() << ","
			<< ci.alpha() << "," << ci.beta() << "," << ci.gamma() << " : " << ci.spacegroup() << "}";
		return os;
	}

	Real A() const { return A_; }
	void A(core::Real Ain) { A_=Ain; }
	Real B() const { return B_; }
	void B(core::Real Bin) { B_=Bin; }
	Real C() const { return C_; }
	void C(core::Real Cin) { C_=Cin; }

	Real alpha() const { return alpha_; }
	void alpha(core::Real alphain) { alpha_=alphain; }
	Real beta() const { return beta_; }
	void beta(core::Real betain) { beta_=betain; }
	Real gamma() const { return gamma_; }
	void gamma(core::Real gammain) { gamma_=gammain; }

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

} // namespace pose
} // namespace core

#endif
