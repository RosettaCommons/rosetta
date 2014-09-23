// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file src/protocols/dna/DNAParameters.hh
/// @brief A class to query for base-paired partners as well as base pair and base step parameters
/// @author Jim Havranek

#ifndef INCLUDED_protocols_dna_DNAParameters_fwd
#define INCLUDED_protocols_dna_DNAParameters_fwd

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace dna {

class DNAParameters;

typedef utility::pointer::shared_ptr< DNAParameters > DNAParametersOP;
typedef utility::pointer::shared_ptr< DNAParameters const > DNAParametersCOP;

} // namespace dna
} // namespace protocols

#endif
