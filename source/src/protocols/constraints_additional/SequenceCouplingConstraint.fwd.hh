// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraints_additional/SequenceCouplingConstraint.fwd.hh
/// @brief  This is a constraint that refers to a core::sequence::SequenceCoupling? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design).
/// @author Hetu Kamisetty

#ifndef INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_fwd_hh
#define INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace constraints_additional {

/// @brief
class SequenceCouplingConstraint;
typedef utility::pointer::shared_ptr< SequenceCouplingConstraint > SequenceCouplingConstraintOP;
typedef utility::pointer::shared_ptr< SequenceCouplingConstraint const > SequenceCouplingConstraintCOP;

} // namespace constraints_additional
} // namespace protocols

#endif // INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_fwd_hh
