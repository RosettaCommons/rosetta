// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/Histogram.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_Histogram_FWD_HH
#define INCLUDED_protocols_recces_Histogram_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace recces {

class Histogram;
typedef utility::pointer::shared_ptr< Histogram > HistogramOP;
typedef utility::pointer::shared_ptr< Histogram const > HistogramCOP;

} //recces
} //protocols

#endif
