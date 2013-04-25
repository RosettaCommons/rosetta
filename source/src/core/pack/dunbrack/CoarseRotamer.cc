// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/CoarseRotamer.cc
/// @brief
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/pack/dunbrack/CoarseRotamer.hh>

namespace core {
namespace pack {
namespace dunbrack {

CoarseRotamer::CoarseRotamer(
	Real prob,
	Size nchi_aa,
	RotVector const & rot_in,
	ChiVector const & chi_mean_in,
	ChiVector const & chi_sdev_in,
	AngleVector const & phi_mean,
	AngleVector const & phi_sdev
) :
	probability_( prob ),
	nchi_aa( nchi_aa ),
	rot_( rot_in ),
	chi_mean_( chi_mean_in ),
	chi_sd_( chi_sdev_in ),
	phi_mean_( phi_mean ),
	phi_sdev_( phi_sdev )
{}

CoarseRotamerSet::CoarseRotamerSet()
{}

void
CoarseRotamerSet::push_back( CoarseRotamerOP rot )
{
	rotamers_.push_back( rot );
}




} // namespace dunbrack
} // namespace scoring
} // namespace core

