// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_closure/ccd/CCDLoopClosureMoverCreator.hh
/// @brief  This class will create upcasted instances of CCDLoopClosureMover for the protocols::moves::MoverFactory
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_closure_ccd_ccd_loop_closure_mover_creator_HH
#define INCLUDED_protocols_loops_loops_closure_ccd_ccd_loop_closure_mover_creator_HH

// unit headers
#include <protocols/moves/MoverCreator.hh>

//C++ headers
#include <string>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

/// @brief %CCDLoopClosureMoverCreator allows the MoverFactory to create a CCDLoopClosureMover instance.
class CCDLoopClosureMoverCreator : public protocols::moves::MoverCreator {

public:
	/// @brief Return a up-casted owning pointer (MoverOP) to the mover.
	virtual protocols::moves::MoverOP create_mover() const;

	/// @brief Return the string identifier for the associated Mover (CCDLoopClosureMover).
	virtual std::string keyname() const;

	/// @brief Static method that returns the keyname for performance reasons.
	static std::string mover_name();
};

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loops_closure_ccd_ccd_loop_closure_mover_creator_HH
