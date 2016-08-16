// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/floppy_tail/FloppyTail_publication.hh
/// @brief FloppyTail extra functions from original publication - this calculates some statistics used for the first published use of the code (Kleiger G, Saha A, Lewis S, Kuhlman B, Deshaies RJ. Rapid E2-E3 assembly and disassembly enable processive ubiquitylation of cullin-RING ubiquitin ligase substrates. Cell. 2009 Nov 25;139(5):957-68. PubMed PMID: 19945379.)
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_floppy_tail_FloppyTail_publication_hh
#define INCLUDED_protocols_floppy_tail_FloppyTail_publication_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace floppy_tail {

/// @details This function is specific to the original system for which this code was written
void create_extra_output( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP score_fxn );

} //floppy_tail
} //protocols

#endif //INCLUDED_protocols_floppy_tail_FloppyTail_publication_hh
