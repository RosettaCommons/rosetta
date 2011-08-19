// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopRelaxMover.fwd.hh
/// @brief forward declaration for LoopRelaxMover class
/// @author James Thompson

#ifndef INCLUDED_protocols_loops_LoopRelaxMover_fwd_hh
#define INCLUDED_protocols_loops_LoopRelaxMover_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace loops {

class LoopRelaxMover;

typedef  utility::pointer::owning_ptr< LoopRelaxMover >  LoopRelaxMoverOP;
typedef  utility::pointer::owning_ptr< LoopRelaxMover const >  LoopRelaxMoverCOP;



class LoopRelaxThreadingMover;

typedef  utility::pointer::owning_ptr< LoopRelaxThreadingMover >  LoopRelaxThreadingMoverOP;
typedef  utility::pointer::owning_ptr< LoopRelaxThreadingMover const >  LoopRelaxThreadingMoverCOP;


} // namespace loops
} // namespace protocols

#endif
