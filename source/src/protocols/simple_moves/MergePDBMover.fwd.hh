// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/MergePdbMoverCreator.hh
/// @brief This class will allign & combine parts of the pdb.
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_MergePdbMover_fwd_hh
#define INCLUDED_protocols_simple_moves_MergePdbMover_fwd_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace simple_moves {

class MergePDBMover;

typedef utility::pointer::shared_ptr< MergePDBMover > MergePDBMoverOP;
typedef utility::pointer::shared_ptr< MergePDBMover const > MergePDBMoverCOP;

} // simple_moves
} // protocols

#endif // INCLUDED_protocols_simple_moves_MergePdbMoverCreator_fwd_hh

