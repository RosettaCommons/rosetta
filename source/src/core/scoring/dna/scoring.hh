// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_dna_scoring_hh
#define INCLUDED_core_scoring_dna_scoring_hh

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.fwd.hh>
// #include <numeric/xyzVector.hh>
// #include <numeric/xyzMatrix.hh>


namespace core { namespace scoring {
namespace dna {

kinematics::Stub
get_base_stub(
	conformation::Residue const & rsd,
	int const strand
);


/* void
base_pair_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< Real > & params // output
); */

/*
void
new_base_step_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< Real > & params // output
); */


} // namespace dna
}} // scoring core

#endif
