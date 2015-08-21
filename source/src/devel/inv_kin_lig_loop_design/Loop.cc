// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Loop.cc
///
/// @brief
/// @author
#include <devel/inv_kin_lig_loop_design/Loop.hh>
#include <devel/inv_kin_lig_loop_design/std_extra.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <core/chemical/ResidueConnection.hh>

namespace devel {

namespace inv_kin_lig_loop_design {

Loop::Loop(Segment const& segment, map<Residue*,Residue*> const& clones ) : lo(0), hi(0), to(0), jumpno(0) {

	lo = find_or_throw(clones, segment.lo_res )->seqpos();
	hi = find_or_throw(clones, segment.hi_res )->seqpos();
	tag = segment.tag;

	if ( segment.type == Segment::ANCHORED_LOOP ) {
		to = lo + segment.nres_pre;
		jumpno = segment.jumpno;
	}

}

} // LoopDesign

} // devel
