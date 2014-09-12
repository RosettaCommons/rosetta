// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/MultiFingerprintMultifunc.hh
/// @brief  Fingerprint comparison multifunction class
/// @author Ragul Gowthaman


#ifndef INCLUDED_protocols_pockets_MultiFingerprintMultifunc_hh
#define INCLUDED_protocols_pockets_MultiFingerprintMultifunc_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>

#include <protocols/pockets/Fingerprint.fwd.hh>
#include <protocols/pockets/MultiFingerprintMultifunc.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace pockets {

/// @brief Atom tree multifunction class
class MultiFingerprintMultifunc : public core::optimization::Multifunc {

public: // Creation

	// c-tor
	MultiFingerprintMultifunc(
		NonPlaidFingerprint & nfp_in,
		PlaidFingerprint & pfp_in,
		core::Real const & missing_point_weight,
		core::Real const & steric_weight,
		core::Real const & extra_point_weight,
		core::Size const nconformers
	);

	/// @brief Destructor
	inline
	virtual
	~MultiFingerprintMultifunc()
	{}

public: // Methods

	// func
	virtual
	core::Real
	operator ()( core::optimization::Multivec const & vars ) const;

	// dfunc
	virtual
	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const;

	using core::optimization::Multifunc::dump;

	/// @brief Error state reached; dump out current pdb.
	virtual
	void
	dump( core::optimization::Multivec const & vars ) const;

private:

  NonPlaidFingerprint & nfp_;
  PlaidFingerprint & pfp_;
	core::Real missing_pt_;
	core::Real steric_;
	core::Real extra_pt_;
	core::Size nconformers_;

}; // MultiFingerprintMultifunc


} // namespace Pockets
} // namespace protocols


#endif // INCLUDED_protocols_pockets_MultiFingerprintMultifunc_HH
