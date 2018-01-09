// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/MakeJunctionsMover.fwd.hh
/// @brief
/// @detailed
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_MakeJunctionsMover_FWD_HH
#define INCLUDED_protocols_simple_moves_MakeJunctionsMover_FWD_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace simple_moves {

class MakeJunctionsMover;
typedef utility::pointer::shared_ptr< MakeJunctionsMover > MakeJunctionsMoverOP;
typedef utility::pointer::shared_ptr< MakeJunctionsMover const > MakeJunctionsMoverCOP;

} //simple_moves
} //protocols

#endif
