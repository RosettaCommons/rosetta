// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.fwd.hh
/// @brief  Defines owning pointers for SimpleCycpepPredictApplication class.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace cyclic_peptide_predict {

class SimpleCycpepPredictApplication; // fwd declaration
typedef utility::pointer::shared_ptr< SimpleCycpepPredictApplication > SimpleCycpepPredictApplicationOP;
typedef utility::pointer::shared_ptr< SimpleCycpepPredictApplication const > SimpleCycpepPredictApplicationCOP;
typedef utility::vector1<SimpleCycpepPredictApplicationOP> SimpleCycpepPredictApplicationOPs;
typedef utility::vector1<SimpleCycpepPredictApplicationCOP> SimpleCycpepPredictApplicationCOPs;

} // namespace cyclic_peptide_predict
} // namespace protocols

#endif // INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_fwd_hh
