// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RNA_SuitePotential.fwd.hh
/// @brief  RNA_SuitePotential potential class forward delcaration
/// @author Fang-Chieh Chou

#ifndef INCLUDED_core_scoring_rna_RNA_SuitePotential_fwd_hh
#define INCLUDED_core_scoring_rna_RNA_SuitePotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace rna {

class RNA_SuitePotential;

typedef  utility::pointer::shared_ptr< RNA_SuitePotential > RNA_SuitePotentialOP;
typedef  utility::pointer::shared_ptr< RNA_SuitePotential const > RNA_SuitePotentialCOP;

} //rna
} //scoring
} //core

#endif
