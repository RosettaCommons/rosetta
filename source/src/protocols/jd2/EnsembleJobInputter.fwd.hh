// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/jd2/EnsembleJobInputter.fwd.hh
/// @brief A Job Inputter for distributing a job based on a set of input structures that make up an ensemble
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_jd2_EnsembleJobInputter_fwd_hh
#define INCLUDED_protocols_jd2_EnsembleJobInputter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace jd2 {

class EnsembleJobInputter;

typedef utility::pointer::shared_ptr< EnsembleJobInputter > EnsembleJobInputterOP;
typedef utility::pointer::shared_ptr< EnsembleJobInputter const > EnsembleJobInputterCOP;



} //protocols
} //jd2


#endif //INCLUDED_protocols_jd2_EnsembleJobInputter_fwd_hh





