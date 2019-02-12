// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/TCRloopRefine.cc
/// @brief TCR cdr loop modeling/refinement
/// @author Ragul Gowthaman (ragul@umd.edu)

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

// Basic headers
#include <core/types.hh>
#include <basic/Tracer.hh>
// Protocol includes
#include <protocols/tcr/TCRloopRefine.fwd.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <protocols/tcr/TCRloopRefine.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/modeling_util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
// Option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/tcrmodel.OptionKeys.gen.hh>
// Utility Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <regex>

static basic::Tracer TR( "protocols.tcr.TCRloopRefine" );

using namespace core;
using namespace basic::options;

namespace protocols {
namespace tcr {

TCRloopRefine::TCRloopRefine( TCRseqInfo const &tcrseq_info) {
	//ts_info_ = tcrseq_info;
	aseqs_ = tcrseq_info.atcr();
	aseqs_ = tcrseq_info.atcr();
	bseqs_ = tcrseq_info.btcr();
	aposi_ = tcrseq_info.aseq_posi();
	bposi_ = tcrseq_info.bseq_posi();
	init();
	set_default();
}

TCRloopRefine::TCRloopRefine( TCRseqInfo::tcrsegs const &atcrseqs,
	TCRseqInfo::tcrsegs const &btcrseqs,
	TCRseqInfo::tcrposi const &acdrposi,
	TCRseqInfo::tcrposi const &bcdrposi ) {
	aseqs_ = atcrseqs;
	bseqs_ = btcrseqs;
	aposi_ = acdrposi;
	bposi_ = bcdrposi;
	init();
	set_default();
}

TCRloopRefineOP TCRloopRefine::clone() const{
	return utility::pointer::make_shared< TCRloopRefine >( * this);
}

TCRloopRefine::~TCRloopRefine() = default;

void TCRloopRefine::init() {
	outname_prefix_ = basic::options::option[basic::options::OptionKeys::out::prefix];
	scorefxn_ = core::scoring::get_score_function();
	ind_remodel_cdr3a_loop_ = option[OptionKeys::tcrmodel::remodel_loop_cdr3a];
	ind_remodel_cdr3b_loop_ = option[OptionKeys::tcrmodel::remodel_loop_cdr3b];
	remodel_cdr3a_loop_ = option[OptionKeys::tcrmodel::remodel_tcr_cdr3a_loop];
	remodel_cdr3b_loop_ = option[OptionKeys::tcrmodel::remodel_tcr_cdr3b_loop];
	remodel_cdr3_loops_ = option[OptionKeys::tcrmodel::remodel_tcr_cdr3_loops];
	refine_cdr3a_loop_ = option[OptionKeys::tcrmodel::refine_tcr_cdr3a_loop];
	refine_cdr3b_loop_ = option[OptionKeys::tcrmodel::refine_tcr_cdr3b_loop];
	refine_cdr3_loops_ = option[OptionKeys::tcrmodel::refine_tcr_cdr3_loops];
	refine_all_cdr_loops_ = option[OptionKeys::tcrmodel::refine_all_tcr_cdr_loops];
	return;
}

void TCRloopRefine::set_default() {
	if ( remodel_cdr3a_loop_ && remodel_cdr3b_loop_ ) {
		remodel_cdr3_loops_ = true;
	}
	if ( refine_cdr3a_loop_ && refine_cdr3b_loop_ ) {
		refine_cdr3_loops_ = true;
	}
	if ( ind_remodel_cdr3a_loop_ ||
			ind_remodel_cdr3b_loop_ ||
			remodel_cdr3a_loop_ ||
			remodel_cdr3b_loop_ ||
			remodel_cdr3_loops_ ||
			refine_cdr3a_loop_ ||
			refine_cdr3b_loop_ ||
			refine_cdr3_loops_ ||
			refine_all_cdr_loops_
			) {
		set_tcr_loop_model(true);
	} else {
		set_tcr_loop_model(false);
	}
	return;
}

void TCRloopRefine::apply(core::pose::Pose &pose_in) {
	if ( tcr_loop_model_ ) {
		(*scorefxn_)(pose_in);
		TR <<"score before loop model : "<<pose_in.energies().total_energies()[ core::scoring::total_score ]<<std::endl;
		core::pose::Pose f_pose = pose_in;
		//Loop remodeling using IndependentLoopMover (default remodel mover : "perturb_kic")
		if ( ind_remodel_cdr3a_loop_ || ind_remodel_cdr3b_loop_ ) {
			if( ! protocols::antibody::grafting::antibody_grafting_usable() ) {
				utility_exit_with_message("ERROR: Your compiler does not have full support for C++11 regex.");
			}
			std::string inpstringa = aseqs_.truncdomain;
			std::smatch matcha;
			std::regex patrna(aseqs_.cdr3);
			std::regex_search ( inpstringa, matcha, patrna );
			core::Size cdr3a_pose_start = matcha.position(0) + 1;
			core::Size cdr3a_pose_end = matcha.position(0) + aseqs_.cdr3.length();
			std::string inpstringb = bseqs_.truncdomain;
			std::smatch matchb;
			std::regex patrnb(bseqs_.cdr3);
			std::regex_search ( inpstringb, matchb, patrnb );
			core::Size cdr3b_pose_start = aseqs_.truncdomain.length() + matchb.position(0) + 1;
			core::Size cdr3b_pose_end = aseqs_.truncdomain.length() + matchb.position(0) + bseqs_.cdr3.length();
			core::pose::Pose ind_tcr_loopmodel_pose = f_pose;
			if ( ind_remodel_cdr3a_loop_ && ind_remodel_cdr3b_loop_ ) {
				ind_remodel_tcr_cdr_loops( ind_tcr_loopmodel_pose, cdr3a_pose_start, cdr3a_pose_end );
				ind_remodel_tcr_cdr_loops( ind_tcr_loopmodel_pose, cdr3b_pose_start, cdr3b_pose_end );
			} else if ( ind_remodel_cdr3a_loop_ ) {
				ind_remodel_tcr_cdr_loops( ind_tcr_loopmodel_pose, cdr3a_pose_start, cdr3a_pose_end );
			} else if ( ind_remodel_cdr3b_loop_ ) {
				ind_remodel_tcr_cdr_loops( ind_tcr_loopmodel_pose, cdr3b_pose_start, cdr3b_pose_end );
			}
			f_pose = ind_tcr_loopmodel_pose;
		}
		//Loop remodeling (and refinement)
		if ( remodel_cdr3_loops_ || remodel_cdr3a_loop_ || remodel_cdr3b_loop_ ) {
			core::pose::Pose tcr_loopmodel_pose = f_pose;
			protocols::loops::Loop cdr3aloop( protocols::loops::Loop(aposi_.cdr3.begin, aposi_.cdr3.end) );
			protocols::loops::Loop cdr3bloop( protocols::loops::Loop(bposi_.cdr3.begin, bposi_.cdr3.end) );
			if ( remodel_cdr3_loops_ ) {
				minimize_tcr_model( tcr_loopmodel_pose );
				protocols::loops::LoopsOP cdr3loops( new protocols::loops::Loops() );
				cdr3loops->add_loop(cdr3aloop);
				cdr3loops->add_loop(cdr3bloop);
				bool refine = false;
				if ( refine_cdr3_loops_ ) refine = true;
				remodel_tcr_cdr_loops( tcr_loopmodel_pose, cdr3loops, refine );
			} else {
				if ( remodel_cdr3a_loop_ ) {
					minimize_tcr_model( tcr_loopmodel_pose );
					protocols::loops::LoopsOP cdr3loops( new protocols::loops::Loops() );
					cdr3loops->add_loop(cdr3aloop);
					bool refine = false;
					if ( refine_cdr3a_loop_ ) refine = true;
					remodel_tcr_cdr_loops( tcr_loopmodel_pose, cdr3loops, refine );
				}
				if ( remodel_cdr3b_loop_ ) {
					minimize_tcr_model( tcr_loopmodel_pose );
					protocols::loops::LoopsOP cdr3loops( new protocols::loops::Loops() );
					cdr3loops->add_loop(cdr3bloop);
					bool refine = false;
					if ( remodel_cdr3b_loop_ ) refine = true;
					remodel_tcr_cdr_loops( tcr_loopmodel_pose, cdr3loops, refine );
				}
			}
			f_pose = tcr_loopmodel_pose;
		}
		//Loop refinement only
		if ( refine_all_cdr_loops_ ) {
			minimize_tcr_model( f_pose );
			refine_all_cdr_loops( f_pose, aposi_, bposi_ );
		} else if ( refine_cdr3_loops_ && !remodel_cdr3_loops_ ) {
			minimize_tcr_model( f_pose );
			refine_cdr3_loops( f_pose, aposi_.cdr3, bposi_.cdr3 );
		} else {
			if ( refine_cdr3a_loop_ && !remodel_cdr3a_loop_ ) {
				refine_cdr_loop( f_pose, aposi_.cdr3 );
			}
			if ( refine_cdr3b_loop_ && !remodel_cdr3b_loop_ ) {
				refine_cdr_loop( f_pose, bposi_.cdr3 );
			}
		}
		pose_in = f_pose;
		(*scorefxn_)(f_pose);
		TR <<"score after loop model : "<<pose_in.energies().total_energies()[ core::scoring::total_score ]<<std::endl;
	}
	return;
}

} // namespace tcr
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

