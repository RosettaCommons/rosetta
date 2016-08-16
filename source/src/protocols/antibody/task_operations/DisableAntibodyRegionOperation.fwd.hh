// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/task_operations/DisableAntibodyRegionOperation.fwd.hh
/// @brief Task operation to disable regions of an antibody
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_task_operations_DisableAntibodyRegionOperationfwd_hh
#define INCLUDED_protocols_antibody_task_operations_DisableAntibodyRegionOperationfwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace antibody {
namespace task_operations {


// Forward
class DisableAntibodyRegionOperation;

typedef utility::pointer::shared_ptr< DisableAntibodyRegionOperation > DisableAntibodyRegionOperationOP;
typedef utility::pointer::shared_ptr< DisableAntibodyRegionOperation const > DisableAntibodyRegionOperationCOP;


} //task_operations
} //antibody
} //protocols

#endif





