// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file InterModelMotifScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/InterModelMotifScorer.hh>

//Core headers
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/motif/util.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

static THREAD_LOCAL basic::Tracer TR("protocols.sewing.scoring.InterModelMotifScorer");

InterModelMotifScorer::InterModelMotifScorer():
	MotifScorer()
{}

///@details Use the negative normalized motif score
core::Real
InterModelMotifScorer::score(
	AssemblyCOP assembly
) {
	return -1.0 * full_motif_score(assembly);
}

///@details use Will's Motif score to calculate the motif score for interactions between
///a given segment and segments from other models. Divide by total number of segments
core::Real
InterModelMotifScorer::full_motif_score(
	AssemblyCOP assembly
) {
	core::Real score = 0.0;
	utility::vector1<SewSegment> const segments = assembly->segments();
	core::Size counter = 0;
	for ( core::Size i=1; i<=segments.size(); ++i ) {
		for ( core::Size res_i = 1; res_i <= segments[i].residues_.size(); ++res_i ) {

			numeric::xyzTransform<core::Real> stub1 = get_stub(segments, i, res_i);
			char ss1 = segments[i].dssp_;
			char aa1 = res_type_set_->name_map(segments[i].residues_[res_i].residue_type_).name1();

			//******DIRTY TEMPORARY HACK WILL ONLY WORK FOR CONTINUOUS ASSEMBLIES WHERE ALL MODELS ARE 3 SEGMENTS LONG*********//
			for ( core::Size j=i+3; j<=segments.size(); ++j ) {
				for ( core::Size res_j = 1; res_j <= segments[j].residues_.size(); ++res_j ) {
					numeric::xyzTransform<core::Real> stub2 = get_stub(segments, j, res_j);
					char ss2 = segments[j].dssp_;
					char aa2 = res_type_set_->name_map(segments[j].residues_[res_j].residue_type_).name1();

					score += get_score(stub1, ss1, aa1, stub2, ss2, aa2);
					++counter;
				}
			}
		}
	}
	if ( counter == 0 ) { return score; }
	return score / counter;
}

} //scoring namespace
} //sewing namespace
} //protocols namespace
