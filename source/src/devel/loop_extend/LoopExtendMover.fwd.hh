// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/LoopExtend/LoopExtendMover.fwd.hh
/// @brief  LoopExtend Mover forward declarations header
/// @author Daniel J. Mandell


#ifndef INCLUDED_devel_loop_extend_LoopExtendMover_fwd_hh
#define INCLUDED_devel_loop_extend_LoopExtendMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace devel {
namespace loop_extend {

//Forwards and OP typedefs
class LoopExtendMover;
typedef utility::pointer::shared_ptr< LoopExtendMover > LoopExtendMoverOP;


}//LoopExtend
}//devel

#endif //INCLUDED_devel_LoopExtend_LoopExtendMover_FWD_HH
