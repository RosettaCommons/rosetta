// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyPartnerMotifScorer.cc
///
/// @brief
/// @author Frank Teets

//Unit headers
#include <protocols/legacy_sewing/scoring/LegacyPartnerMotifScorer.hh>

//Core headers
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/motif/util.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/Pose.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace legacy_sewing  {
namespace scoring {

static THREAD_LOCAL basic::Tracer TR("protocols.legacy_sewing.scoring.LegacyPartnerMotifScorer");

LegacyPartnerMotifScorer::LegacyPartnerMotifScorer():
	LegacyMotifScorer()
{}

///@details Use the negative normalized motif score
core::Real
LegacyPartnerMotifScorer::score(
	AssemblyCOP assembly
) {
	return -1.0 * interface_motif_score(assembly);
}

///@details use Will's Motif score to calculate the motif score for the entire Assembly
core::Real
LegacyPartnerMotifScorer::interface_motif_score(
	AssemblyCOP assembly
) {
	core::Real score = 0.0;

	utility::vector1<SewSegment> const segments = assembly->segments();
	core::pose::PoseOP const partner = assembly->get_partner();
	if ( !partner ) {
		return score;
	}
	if ( partner->secstruct() == "" ) {
		core::scoring::dssp::Dssp dssp(*partner);
		dssp.insert_ss_into_pose(*partner);
	}

	for ( core::Size i=1; i<=segments.size(); ++i ) {
		for ( core::Size res_i = 1; res_i <= segments[i].residues_.size(); ++res_i ) {

			numeric::xyzTransform<core::Real> stub1 = get_stub(segments, i, res_i);
			char ss1 = segments[i].dssp_;
			char aa1 = res_type_set_->name_map(segments[i].residues_[res_i].residue_type_).name1();

			for ( core::Size j=1; j<=partner->size(); ++j ) {

				numeric::xyzTransform<core::Real> stub2 = core::pose::motif::get_backbone_reference_frame(*partner, j);
				char ss2 = partner->secstruct(j);
				char aa2 = partner->residue(j).name1();

				score += get_score(stub1, ss1, aa1, stub2, ss2, aa2);
			}
		}
	}
	return score / (partner->size() + assembly->total_residue());
}

} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace
