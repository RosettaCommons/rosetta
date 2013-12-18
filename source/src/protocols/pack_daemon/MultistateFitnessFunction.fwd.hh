// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/MultistateFitnessFunction.fwd.hh
/// @brief  forward declaration of class MultistateFitnessFunction
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_MultistateFitnessFunction_fwd_hh
#define INCLUDED_protocols_pack_daemon_MultistateFitnessFunction_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace pack_daemon {

class TopEntitySet;

typedef utility::pointer::owning_ptr< TopEntitySet > TopEntitySetOP;
typedef utility::pointer::owning_ptr< TopEntitySet const > TopEntitySetCOP;


class MultistateFitnessFunction;

typedef utility::pointer::owning_ptr< MultistateFitnessFunction > MultistateFitnessFunctionOP;
typedef utility::pointer::owning_ptr< MultistateFitnessFunction const > MultistateFitnessFunctionCOP;

class MPIMultistateFitnessFunction;

typedef utility::pointer::owning_ptr< MPIMultistateFitnessFunction > MPIMultistateFitnessFunctionOP;
typedef utility::pointer::owning_ptr< MPIMultistateFitnessFunction const > MPIMultistateFitnessFunctionCOP;


}
}

#endif
