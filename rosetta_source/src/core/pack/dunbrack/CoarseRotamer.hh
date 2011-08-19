// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/CoarseRotamer.hh
/// @brief
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_CoarseRotamer_hh
#define INCLUDED_core_pack_dunbrack_CoarseRotamer_hh

// Unit headers
#include <core/pack/dunbrack/CoarseRotamer.fwd.hh>

// Package headers
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <utility/vector1_bool.hh>
#include <iterator>



namespace core {
namespace pack {
namespace dunbrack {

class CoarseRotamer : public utility::pointer::ReferenceCount
{
public:
	CoarseRotamer(
		Real prob,
		Size nchi_aa,
		RotVector const & rot_in,
		ChiVector const & chi_mean_in,
		ChiVector const & chi_sdev_in,
		AngleVector const & phi_mean,
		AngleVector const & phi_sdev
	);

private:
	Real probability_;
	Size nchi_aa;
	RotVector rot_;
	ChiVector chi_mean_;
	ChiVector chi_sd_;
	AngleVector phi_mean_;
	AngleVector phi_sdev_;

};

class CoarseRotamerSet : public utility::pointer::ReferenceCount
{
public:
	CoarseRotamerSet();

	void
	push_back( CoarseRotamerOP rot );

private:
	utility::vector1< CoarseRotamerOP > rotamers_;
};



} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_CoarseSingleResidueLibrary_FWD_HH
