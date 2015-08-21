// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopMover_CCDCreator.hh
/// @brief  Header for LoopMover_CCDCreator
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCDCreator_hh
#define INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCDCreator_hh

// Unit Headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

/// @brief creator for the LoopMover_Refine_CCDCreator class
class LoopMover_Refine_CCDCreator : public moves::MoverCreator
{
public:
	virtual ~LoopMover_Refine_CCDCreator();

	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;

};

} //namespace refine
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCDCreator_hh
