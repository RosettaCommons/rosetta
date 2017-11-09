// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyMotifScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/legacy_sewing/scoring/LegacyMotifScorer.hh>

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
namespace legacy_sewing  {
namespace scoring {

static THREAD_LOCAL basic::Tracer TR("protocols.legacy_sewing.scoring.LegacyMotifScorer");

LegacyMotifScorer::LegacyMotifScorer():
	mman_(*core::scoring::motif::MotifHashManager::get_instance()),
	res_type_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) )
{}

///@details Use the negative normalized motif score
core::Real
LegacyMotifScorer::score(
	AssemblyCOP assembly
) {
	return -1.0 * norm_motif_score(assembly);
}


///@details Motif score of the entire Assembly divided by total residue
core::Real
LegacyMotifScorer::norm_motif_score(
	AssemblyCOP assembly
) {
	return full_motif_score(assembly)/assembly->total_residue();
}


///@details use Will's Motif score to calculate the motif score for the entire Assembly
core::Real
LegacyMotifScorer::full_motif_score(
	AssemblyCOP assembly
) {
	core::Real score = 0.0;
	utility::vector1<SewSegment> const segments = assembly->segments();
	for ( core::Size i=1; i<=segments.size(); ++i ) {
		for ( core::Size res_i = 1; res_i <= segments[i].residues_.size(); ++res_i ) {

			numeric::xyzTransform<core::Real> stub1 = get_stub(segments, i, res_i);
			char ss1 = segments[i].dssp_;
			char aa1 = res_type_set_->name_map(segments[i].residues_[res_i].residue_type_).name1();

			for ( core::Size j=i; j<=segments.size(); ++j ) {
				for ( core::Size res_j = 1; res_j <= segments[j].residues_.size(); ++res_j ) {
					if ( i == j && res_i == res_j ) { continue; }

					numeric::xyzTransform<core::Real> stub2 = get_stub(segments, j, res_j);
					char ss2 = segments[j].dssp_;
					char aa2 = res_type_set_->name_map(segments[j].residues_[res_j].residue_type_).name1();

					score += get_score(stub1, ss1, aa1, stub2, ss2, aa2);
				}
			}
		}
	}
	return score;
}

numeric::xyzTransform<core::Real>
LegacyMotifScorer::get_stub(
	utility::vector1<SewSegment> const & segments,
	core::Size segment_num,
	core::Size resnum
) const {
	return core::pose::motif::get_backbone_reference_frame(
		segments[segment_num].residues_[resnum].basis_atoms_[1].coords_,
		segments[segment_num].residues_[resnum].basis_atoms_[2].coords_,
		segments[segment_num].residues_[resnum].basis_atoms_[3].coords_);
}

core::Real
LegacyMotifScorer::get_score(
	numeric::xyzTransform<core::Real> stub1,
	char ss1,
	char aa1,
	numeric::xyzTransform<core::Real> stub2,
	char ss2,
	char aa2
) const {
	core::Real score = 0;
	core::scoring::motif::XformScoreCOP xs_bb_fxn1 = mman_.get_xform_score_BB_BB(ss1,ss2,aa1,aa2);
	if ( !xs_bb_fxn1 ) {
		std::stringstream err;
		err << "Null XformScore pointer!" << std::endl;
		err << "ss1 " << ss1 << " aa1" << aa1 << std::endl;
		err << "ss2 " << ss2 << " aa2" << aa2 << std::endl;
		utility_exit_with_message(err.str());
	}
	numeric::xyzTransform<core::Real> const Xbb = stub1.inverse() * stub2;
	if ( Xbb.x() < 16 && Xbb.x() > -16
			&& Xbb.y() < 16 && Xbb.y() > -16
			&& Xbb.z() < 16 && Xbb.z() > -16
			) {
		score += xs_bb_fxn1->score_of_bin(Xbb);
	}

	core::scoring::motif::XformScoreCOP xs_bb_fxn2 = mman_.get_xform_score_BB_BB(ss2,ss1,aa2,aa1);
	if ( !xs_bb_fxn2 ) {
		std::stringstream err;
		err << "Null XformScore pointer!" << std::endl;
		err << "ss1 " << ss1 << " aa1" << aa1 << std::endl;
		err << "ss2 " << ss2 << " aa2" << aa2 << std::endl;
		utility_exit_with_message(err.str());
	}
	numeric::xyzTransform<core::Real> const Xbbi = Xbb.inverse();
	if ( Xbbi.x() < 16 && Xbbi.x() > -16
			&& Xbbi.y() < 16 && Xbbi.y() > -16
			&& Xbbi.z() < 16 && Xbbi.z() > -16
			) {
		score += xs_bb_fxn2->score_of_bin(Xbbi);
	}
	return score;
}

void
LegacyMotifScorer::dump_motif(
	AssemblyCOP assembly
) const {
	core::pose::Pose allpose(assembly->to_pose(core::chemical::FA_STANDARD, false));
	core::pose::Pose assembledpose = *(allpose.split_by_chain(1));
	core::scoring::dssp::Dssp dssp(assembledpose);
	dssp.insert_ss_into_pose(assembledpose);
	core::scoring::motif::ResPairMotifQuery opt(assembledpose);
	opt.interface_only() = false;
	core::scoring::motif::MotifHits hits;
	core::scoring::motif::MotifHashManager::get_instance()->get_matching_motifs(opt,hits);
	if ( hits.size() ) {
		std::string outfile = "motif_test.pdb";
		std::cout << "dump " << hits.size() << " (before pruning) hits to " << outfile << std::endl;
		utility::io::ozstream pdbout( outfile );
		pdbout << "MODEL MAIN" << std::endl;
		assembledpose.dump_pdb(pdbout);
		pdbout << "ENDMDL" << std::endl;
		hits.dump_motifs_pdb(pdbout);
		pdbout.close();
	}
}

} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace
