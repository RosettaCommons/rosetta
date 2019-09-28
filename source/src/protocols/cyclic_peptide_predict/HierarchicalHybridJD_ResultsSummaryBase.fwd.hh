// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_ResultsSummaryBase.fwd.hh
/// @brief A pure virtual base class for the helper classes used by the HierarchicalHybridJDApplication
/// class for transmitting information up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_ResultsSummaryBase_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_ResultsSummaryBase_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace cyclic_peptide_predict {

class HierarchicalHybridJD_ResultsSummaryBase;

using HierarchicalHybridJD_ResultsSummaryBaseOP = utility::pointer::shared_ptr< HierarchicalHybridJD_ResultsSummaryBase >;
using HierarchicalHybridJD_ResultsSummaryBaseCOP = utility::pointer::shared_ptr< HierarchicalHybridJD_ResultsSummaryBase const >;

} //cyclic_peptide_predict
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_ResultsSummaryBase_fwd_hh
