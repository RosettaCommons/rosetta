// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/comparative_modeling/coord_util.hh
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_coord_util_hh
#define INCLUDED_protocols_comparative_modeling_coord_util_hh

#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>


namespace protocols {
namespace comparative_modeling {

void gather_coords(
	core::pose::Pose const & model,
	core::pose::Pose const & native,
	core::sequence::SequenceAlignment const & aln,
	int & natoms,
	ObjexxFCL::FArray2D< core::Real > & p1a,
	ObjexxFCL::FArray2D< core::Real > & p2a,
	std::string const & atom_name = "CA"
);

} // comparative_modeling
} // protocols

#endif
