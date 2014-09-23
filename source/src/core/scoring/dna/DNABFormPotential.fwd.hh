// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/dna/DNABFormPotential.fwd.hh
/// @brief  DNA B-form specific torsion potential class forward delcaration
/// @author Jim Havranek

#ifndef INCLUDED_core_scoring_dna_DNABFormPotential_FWD_HH
#define INCLUDED_core_scoring_dna_DNABFormPotential_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace scoring {
namespace dna {

class DNABFormPotential;
typedef utility::pointer::shared_ptr< DNABFormPotential > DNABFormPotentialOP;
typedef utility::pointer::shared_ptr< DNABFormPotential const > DNABFormPotentialCOP;

class TorsionFourierComponent;
typedef utility::pointer::shared_ptr< TorsionFourierComponent > TorsionFourierComponentOP;
typedef utility::pointer::shared_ptr< TorsionFourierComponent const > TorsionFourierComponentCOP;

}
}
}

#endif
