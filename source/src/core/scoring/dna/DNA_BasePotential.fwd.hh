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

#ifndef INCLUDED_core_scoring_dna_DNA_BasePotential_fwd_hh
#define INCLUDED_core_scoring_dna_DNA_BasePotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace dna {

class DNA_BasePotential;

typedef  utility::pointer::shared_ptr< DNA_BasePotential > DNA_BasePotentialOP;
typedef  utility::pointer::shared_ptr< DNA_BasePotential const > DNA_BasePotentialCOP;

} // namespace dna
} // scoring
} // core

#endif
