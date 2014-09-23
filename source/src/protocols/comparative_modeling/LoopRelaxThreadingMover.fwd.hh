// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/LoopRelaxThreadingMover.fwd.hh
/// @brief forward declaration for LoopRelaxThreadingMover class
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_LoopRelaxThreadingMover_fwd_hh
#define INCLUDED_protocols_comparative_modeling_LoopRelaxThreadingMover_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace comparative_modeling {

class LoopRelaxThreadingMover;
typedef utility::pointer::shared_ptr< LoopRelaxThreadingMover > LoopRelaxThreadingMoverOP;
typedef utility::pointer::shared_ptr< LoopRelaxThreadingMover const > LoopRelaxThreadingMoverCOP;

} // namespace comparative_modeling
} // namespace protocols

#endif
