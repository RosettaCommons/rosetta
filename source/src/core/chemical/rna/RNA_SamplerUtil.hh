// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rna/RNA_SamplerUtil.hh
/// @brief
/// @author Rhiju

#ifndef INCLUDED_core_chemical_rna_RNA_SamplerUtil_hh
#define INCLUDED_core_chemical_rna_RNA_SamplerUtil_hh

#include <core/types.hh>
#include <utility/vector1.hh>

namespace core {
namespace chemical {
namespace rna {

void add_values_from_center(
	utility::vector1<core::Real> & torsions,
	Real const center,
	Real const max_range,
	Real const bin_size );

utility::vector1< Real >
get_full_torsions( Real const bin_size = 20.0 );

utility::vector1< Real >
get_epsilon_torsions( Real const delta,
	bool const extra_epsilon,
	Real const bin_size /* = 20.0*/ );

utility::vector1< Real >
get_epsilon_torsions( bool const north_pucker,
	bool const extra_epsilon = true,
	Real const bin_size = 20.0 );


} //ns rna
} //ns chemical
} //ns core

#endif
