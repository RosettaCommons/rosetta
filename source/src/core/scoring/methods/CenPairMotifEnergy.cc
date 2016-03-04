// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/motif/MotifScore.cc
/// @brief  Will's motif score to determine how well packed the protein core is
/// @author TJ Brunette
///
// Unit headers
#include <core/scoring/methods/CenPairMotifEnergyCreator.hh>
#include <core/scoring/methods/CenPairMotifEnergy.hh>


// Project Headers

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/EnergyMap.hh>

#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/motif/reference_frames.hh>
#include <numeric/xyzTransform.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <string>

namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.CenPairMotifEnergy" );

using namespace core::scoring::motif;

methods::EnergyMethodOP
CenPairMotifEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new CenPairMotifEnergy );
}

ScoreTypes
CenPairMotifEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cen_pair_motifs );
	return sts;
}

CenPairMotifEnergy::CenPairMotifEnergy():
	parent(methods::EnergyMethodCreatorOP( new CenPairMotifEnergyCreator ) )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	mman_ = core::scoring::motif::MotifHashManager::get_instance();
	if ( option[OptionKeys::score::motif_residues].user() ) {
		aalist_ =option[OptionKeys::score::motif_residues]();
	}
}

//I have implemented this as a WholeStructureEnergy because the DSSP calls would waste time. But it may be useful in the future to develop this term over each residue
void CenPairMotifEnergy::finalize_total_energy( pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::kinematics::FoldTree;
	double score = 0.0;
	Real max_motif_per_res = option[OptionKeys::score::max_motif_per_res]();
	Size nres1 = pose.n_residue(), nres2=nres1;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		nres2 = core::pose::symmetry::symmetry_info(pose)->num_total_residues_without_pseudo();
		if ( option[OptionKeys::score::motif_ignore_symmmetry]() ) {
			nres2= core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		}
	}
	if ( pose.total_residue() == 0 ) {
		totals[cen_pair_motifs] = 0;
	} else {
		const FoldTree& tree = pose.fold_tree();
		core::scoring::dssp::Dssp dssp( pose );
		dssp.dssp_reduced();
		if ( aalist_.size()>0 ) {
			for ( size_t ii = 1; ii <= aalist_.size(); ++ii ) {
				Size ir = aalist_[ii];
				if ( ir >  pose.total_residue() ) {
					TR << "Warning" << ir << "is above total residue" << pose.total_residue() << std::endl;
				} else {
					Xform const ibb_stub = core::pose::motif::get_backbone_reference_frame(pose,ir);
					char ss1 = dssp.get_dssp_secstruct( ir );
					char aa1 = pose.residue(ir).name1();
					for ( size_t jr = 1; jr <= nres2; ++jr ) {// was ir+1 so it only scanned residues once. but this will scan it twice.
						if ( ir != jr && !tree.is_jump_point(ir) && !tree.is_jump_point(jr) ) {
							Real dist = pose.residue(ir).xyz("CA").distance(pose.residue(jr).xyz("CA"));
							if ( dist < 12 ) {
								char ss2 = dssp.get_dssp_secstruct( jr );
								char aa2 = pose.residue(jr).name1();
								//std::cout << ss1 << ss2 << " " << aa1 << aa2 << "dist" << dist <<  std::endl;
								Xform const jbb_stub = core::pose::motif::get_backbone_reference_frame(pose,jr);
								Xform const Xbb = ibb_stub.inverse() * jbb_stub;
								core::scoring::motif::XformScoreCOP xs_bb_fxn1(mman_->get_xform_score_BB_BB(ss1,ss2,aa1,aa2));
								core::scoring::motif::XformScoreCOP xs_bb_fxn2(mman_->get_xform_score_BB_BB(ss2,ss1,aa2,aa1));
								Real tmpScore = 0;
								if ( xs_bb_fxn1 != NULL ) {
									tmpScore += xs_bb_fxn1->score_of_bin(Xbb);
									tmpScore += xs_bb_fxn2->score_of_bin(Xbb.inverse());
								}
								if ( tmpScore>max_motif_per_res ) {
									tmpScore = max_motif_per_res;
								}
								score+=tmpScore;
							}
						}
					}
				}
			}
			totals[cen_pair_motifs] =(-1*score)/(Real)(aalist_.size());
		} else {
			for ( size_t ir = 1; ir <= nres1; ++ir ) {
				Xform const ibb_stub = core::pose::motif::get_backbone_reference_frame(pose,ir);
				char ss1 = dssp.get_dssp_secstruct( ir );
				char aa1 = pose.residue(ir).name1();
				for ( size_t jr = 1; jr <= nres2; ++jr ) {
					if ( ir != jr && !tree.is_jump_point(ir) && !tree.is_jump_point(jr) ) {
						Real dist = pose.residue(ir).xyz("CA").distance(pose.residue(jr).xyz("CA"));
						if ( dist < 12 ) {
							char ss2 = dssp.get_dssp_secstruct( jr );
							char aa2 = pose.residue(jr).name1();
							//std::cout << ss1 << ss2 << " " << aa1 << aa2 << "dist" << dist <<  std::endl;
							Xform const jbb_stub = core::pose::motif::get_backbone_reference_frame(pose,jr);
							Xform const Xbb = ibb_stub.inverse() * jbb_stub;
							core::scoring::motif::XformScoreCOP xs_bb_fxn1(mman_->get_xform_score_BB_BB(ss1,ss2,aa1,aa2));
							core::scoring::motif::XformScoreCOP xs_bb_fxn2(mman_->get_xform_score_BB_BB(ss2,ss1,aa2,aa1));
							Real tmpScore = 0;
							if ( xs_bb_fxn1 != NULL ) {
								tmpScore += xs_bb_fxn1->score_of_bin(Xbb);
								tmpScore += xs_bb_fxn2->score_of_bin(Xbb.inverse());
							}
							if ( tmpScore>max_motif_per_res ) {
								tmpScore = max_motif_per_res;
							}
							score+=tmpScore;
						}
					}
				}
			}
			totals[cen_pair_motifs] =(-1*score)/(Real)nres1;
		}
	}
}


core::Size
CenPairMotifEnergy::version() const
{
	return 1; // Initial versioning
}

}//motif
}//scoring
}//core
