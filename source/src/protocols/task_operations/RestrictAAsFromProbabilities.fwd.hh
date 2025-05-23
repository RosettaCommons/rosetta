// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/RestrictAAsFromProbabilities.fwd.hh
/// @brief A class to restrict designable amino acids depending on the probabilities from a PerResidueProbabilitiesMetric
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifndef INCLUDED_protocols_task_operations_RestrictAAsFromProbabilities_fwd_hh
#define INCLUDED_protocols_task_operations_RestrictAAsFromProbabilities_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace task_operations {

class RestrictAAsFromProbabilities;

using RestrictAAsFromProbabilitiesOP = utility::pointer::shared_ptr< RestrictAAsFromProbabilities >;
using RestrictAAsFromProbabilitiesCOP = utility::pointer::shared_ptr< RestrictAAsFromProbabilities const >;

} //task_operations
} //protocols

#endif //INCLUDED_protocols_task_operations_RestrictAAsFromProbabilities_fwd_hh
