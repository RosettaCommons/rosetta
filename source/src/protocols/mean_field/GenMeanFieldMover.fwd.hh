// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/mean_field/GenMeanFieldMover.fwd.hh
/// @brief  forward declaration for GenMeanFieldMover
/// @author  Aliza Rubenstein aliza.rubenstein@gmail.com


#ifndef INCLUDED_protocols_mean_field_GenMeanFieldMover_fwd_hh
#define INCLUDED_protocols_mean_field_GenMeanFieldMover_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace mean_field {

// Forward
class GenMeanFieldMover;

// Types
typedef utility::pointer::shared_ptr< GenMeanFieldMover >  GenMeanFieldMoverOP;
typedef utility::pointer::shared_ptr< GenMeanFieldMover const >  GenMeanFieldMoverCOP;

} // namespace mean_field
} // namespace protocols

#endif
