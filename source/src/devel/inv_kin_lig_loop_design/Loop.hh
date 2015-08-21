// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Loop.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_LOOP_HH
#define DEVEL_INVKINLIGLOOPDESIGN_LOOP_HH

#include <vector>
#include <map>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace devel {

namespace inv_kin_lig_loop_design {

using namespace std;
using core::conformation::Residue;
using utility::tag::TagCOP;

// =================================================
// ==================== Segment ====================
// =================================================

struct Segment {
	TagCOP tag; // the tag which this segment represents
	core::pose::PoseOP pose; // the pose containing this segment
	core::conformation::Residue *lo_res, *hi_res, *from_res, *to_res; // the beginning and ending residues of the indel, inclusive, NB: to signal a one residue deletion, set begin=end
	vector< core::chemical::AA > aas;

	enum { ORIGINAL = 1, LOOP, ANCHORED_LOOP, DELETION, LIGAND }; //
	int type; // the type of this Segment
	int nres_pre,nres_post; // LOOP: number of residues before and after the chainbreak, ANCHORED_LOOP: number of residues between [lo,to) and (to,hi]
	int jumpno; // ANCHORED_LOOP: the jump assigned to the jump between from & to

	Segment() : lo_res(0), hi_res(0), from_res(0), to_res(0), type(0), /*, lo(0), hi(0),*/nres_pre(0), nres_post(0), jumpno(0) {}

	int size() { return hi_res->seqpos() - lo_res->seqpos() + 1; }

}; // class Indel

// ==============================================
// ==================== Loop ====================
// ==============================================

struct Loop { // NB: this has to be a different data structure than Segment, since pointers may change during design
	int lo,hi,to,jumpno;
	TagCOP tag;
	Loop() : lo(0), hi(0), to(0), jumpno(0) {}
	Loop( Segment const& segment, map<Residue*,Residue*> const& clones );
};

} // LoopDesign

} // devel

#endif // DEVEL_LOOPDESIGN_LOOP_HH
