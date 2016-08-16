// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/comparative_modeling/GenericJobInputter.fwd.hh
/// @brief  header file for ThreadingJobInputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Oliver Lange

#ifndef INCLUDED_protocols_comparative_modeling_ThreadingJobInputter_fwd_hh
#define INCLUDED_protocols_comparative_modeling_ThreadingJobInputter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace comparative_modeling {

class ThreadingJobInputter;
typedef utility::pointer::shared_ptr< ThreadingJobInputter > ThreadingJobInputterOP;
typedef utility::pointer::shared_ptr< ThreadingJobInputter const > ThreadingJobInputterCOP;

}//comparative_modeling
}//protocols

#endif //INCLUDED_protocols_comparative_modeling_ThreadingJobInputter_FWD_HH
