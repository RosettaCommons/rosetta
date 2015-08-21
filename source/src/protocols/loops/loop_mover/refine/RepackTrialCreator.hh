// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/RepackTrialCreator.hh
/// @brief  Header for LoopMover_KICCreator
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )

#ifndef INCLUDED_protocols_loops_loop_mover_refine_RepackTrialCreator_HH
#define INCLUDED_protocols_loops_loop_mover_refine_RepackTrialCreator_HH

// Unit Headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

/// @brief creator for the RepackTrial class
class RepackTrialCreator : public moves::MoverCreator {
public:
	RepackTrialCreator() {}
	virtual ~RepackTrialCreator();

	virtual moves::MoverOP create_mover() const;

	virtual std::string keyname() const;

};

} //namespace refine
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif // INCLUDED_protocols_loops_loop_mover_refine_RepackTrialCreator_HH
