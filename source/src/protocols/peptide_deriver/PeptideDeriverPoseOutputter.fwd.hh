// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverPoseOutputter.fwd.hh
/// @brief outputs poses at different points of the Peptiderive protocol, according to the given set of options
/// @author orlypolo (orlymarcu@gmail.com)

#ifndef INCLUDED_protocols_peptide_deriver_PeptideDeriverPoseOutputter_fwd_hh
#define INCLUDED_protocols_peptide_deriver_PeptideDeriverPoseOutputter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace peptide_deriver {

class PeptideDeriverPoseOutputter;

typedef utility::pointer::shared_ptr< PeptideDeriverPoseOutputter > PeptideDeriverPoseOutputterOP;
typedef utility::pointer::shared_ptr< PeptideDeriverPoseOutputter const > PeptideDeriverPoseOutputterCOP;

} //protocols
} //peptide_deriver

#endif //INCLUDED_protocols_peptide_deriver_PeptideDeriverPoseOutputter_fwd_hh
