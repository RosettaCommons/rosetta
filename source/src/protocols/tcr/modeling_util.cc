// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/modeling_util.cc
/// @brief modeling utility functions for protocols/tcr
/// @author Ragul (ragul@umd.edu)


#include <protocols/tcr/modeling_util.hh>
#include <protocols/tcr/util.hh>
#include <protocols/tcr/grafting_util.hh>
#include <protocols/tcr/TCRloopRefine.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRmodel.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/antibody/util.hh>
#include <protocols/grafting/util.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/types.hh>
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
//option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.tcr.modeling_util" );

namespace protocols {
namespace tcr {

void match_template_and_target_sequence( core::pose::Pose const &scafold, std::string const &targetseq, core::pose::Pose &posecopy ) {
	std::string templateseq = scafold.sequence();
	if ( templateseq.length() == targetseq.length() ) {
		for ( core::Size i=0; i<targetseq.size(); ) {
			char const aa( targetseq[i] );
			core::Size j = ++i;//Pose numbering starts at 1
			if ( j == targetseq.size() ) {
				//protocols::simple_moves::MutateResidue mutate( j, aa );
				std::string three_letter_aa = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( aa ) );
				std::string aatermi = three_letter_aa + ":CtermProteinFull";
				protocols::simple_moves::MutateResidue mutate( j, aatermi );
				mutate.set_update_polymer_dependent( true );
				mutate.apply( posecopy );
			} else {
				if ( scafold.residue( j ).name1() != aa ) {
					protocols::simple_moves::MutateResidue mutate( j, aa );
					mutate.apply( posecopy );
				}
			}
		}
	} else {
		TR <<templateseq.length() <<" "<< targetseq.length()<<std::endl;
		TR <<templateseq <<" "<< targetseq<<std::endl;
		utility_exit_with_message("Error template and target sequence lengths are not same.");
	}
	return;
}


void orient_tcr_chain( core::pose::Pose &fr_piece, core::pose::Pose const &ori_tmplt_chain, TCRseqInfo::tcrposi const &selepos) {
	utility::vector1< std::pair<core::Size, core::Size> > oripos = get_cdr_pdb_pos_from_input_pose(ori_tmplt_chain, selepos);
	utility::vector1< std::pair<core::Size, core::Size> > frpos = get_cdr_pdb_pos_from_input_pose(fr_piece, selepos);
	core::id::AtomID_Map< core::id::AtomID > ori_atom_map;
	core::pose::initialize_atomid_map( ori_atom_map, fr_piece, core::id::AtomID::BOGUS_ATOM_ID() );
	for ( core::Size i = 1; i <= oripos.size(); ++i ) {
		if ( ((oripos[i].first - oripos[i].second) - (frpos[i].first - frpos[i].second)) < 0.001 ) {
			for ( core::Size j=oripos[i].first,k=frpos[i].first; ( j<=oripos[i].second && k<=frpos[i].second); j++,k++ ) {
				core::id::AtomID const caid1( fr_piece.residue(i).atom_index("CA"), k );
				core::id::AtomID const caid2( ori_tmplt_chain.residue(i).atom_index("CA"), j );
				ori_atom_map.set( caid1, caid2 );
				core::id::AtomID const nid1( fr_piece.residue(i).atom_index("N"), k );
				core::id::AtomID const nid2( ori_tmplt_chain.residue(i).atom_index("N"), j );
				ori_atom_map.set( nid1, nid2 );
				core::id::AtomID const cid1( fr_piece.residue(i).atom_index("C"), k );
				core::id::AtomID const cid2( ori_tmplt_chain.residue(i).atom_index("C"), j );
				ori_atom_map.set( cid1, cid2 );
			}
		} else {
			TR <<"Orientation template segment does not match with framework template"<<std::endl;
			utility_exit_with_message("Orientation template segment does not match with framework template");
		}
	}
	core::scoring::superimpose_pose( fr_piece, ori_tmplt_chain, ori_atom_map );
	return;
}

void refine_tcr_cdr_loops ( core::pose::Pose &tcr_pose, protocols::loops::LoopsOP const &cdrloops ) {
	cdrloops->auto_choose_cutpoints( tcr_pose );
	protocols::loops::loop_mover::refine::LoopMover_Refine_KIC loopmover4 (protocols::loops::loop_mover::refine::LoopMover_Refine_KIC(cdrloops, core::scoring::get_score_function()));
	loopmover4.apply( tcr_pose );
	return;
}

void refine_all_cdr_loops ( core::pose::Pose &tcr_pose, TCRseqInfo::tcrposi const &a_posi, TCRseqInfo::tcrposi const &b_posi ) {
	protocols::loops::LoopsOP cdrloops( new protocols::loops::Loops );
	protocols::loops::Loop cdr1aloop( protocols::loops::Loop(a_posi.cdr1.begin, a_posi.cdr1.end) );
	protocols::loops::Loop cdr3aloop( protocols::loops::Loop(a_posi.cdr3.begin, a_posi.cdr3.end) );
	protocols::loops::Loop cdr1bloop( protocols::loops::Loop(b_posi.cdr1.begin, b_posi.cdr1.end) );
	protocols::loops::Loop cdr3bloop( protocols::loops::Loop(b_posi.cdr3.begin, b_posi.cdr3.end) );
	cdrloops->add_loop(cdr1aloop);
	cdrloops->add_loop(cdr3aloop);
	cdrloops->add_loop(cdr1bloop);
	cdrloops->add_loop(cdr3bloop);
	refine_tcr_cdr_loops ( tcr_pose, cdrloops );
	return;
}

void refine_cdr3_loops ( core::pose::Pose &tcr_pose, antibody::grafting::CDR_Bounds const &acdr3, antibody::grafting::CDR_Bounds const &bcdr3 ) {
	protocols::loops::LoopsOP cdrloops( new protocols::loops::Loops );
	protocols::loops::Loop cdr3aloop( protocols::loops::Loop(acdr3.begin, acdr3.end) );
	protocols::loops::Loop cdr3bloop( protocols::loops::Loop(bcdr3.begin, bcdr3.end) );
	cdrloops->add_loop(cdr3aloop);
	cdrloops->add_loop(cdr3bloop);
	refine_tcr_cdr_loops ( tcr_pose, cdrloops );
	return;
}

void refine_cdr_loop ( core::pose::Pose &tcr_pose, antibody::grafting::CDR_Bounds const &cdr) {
	protocols::loops::LoopsOP cdrloops( new protocols::loops::Loops );
	protocols::loops::Loop cdrloop( protocols::loops::Loop(cdr.begin, cdr.end) );
	cdrloops->add_loop(cdrloop);
	refine_tcr_cdr_loops ( tcr_pose, cdrloops );
	return;
}

void remodel_tcr_cdr_loops ( core::pose::Pose &tcr_pose, protocols::loops::LoopsOP const &cdr3loops, bool refine ) {
	//define cdr3 loops
	using namespace basic::options;
	if ( option[ basic::options::OptionKeys::loops::extended ]() ) {
		cdr3loops->set_extended( true );
	} else {
		cdr3loops->set_extended( false );
	}
	cdr3loops->auto_choose_cutpoints( tcr_pose );
	//utility::vector1<bool> allow_chi_copy( tcr_pose.total_residue(), true );
	//for ( core::Size i=bpose_posi.cdr3start; i<=bpose_posi.cdr3end; i++ ) {
	// allow_chi_copy[i] = false;
	//}
	//protocols::simple_moves::ReturnSidechainMover recover_sidechains( tcr_pose, allow_chi_copy );
	core::kinematics::FoldTree original_ftree( tcr_pose.fold_tree() );
	core::kinematics::FoldTree new_ftree;
	protocols::loops::fold_tree_from_loops( tcr_pose, *cdr3loops, new_ftree, false );
	tcr_pose.fold_tree( new_ftree );
	protocols::loops::add_cutpoint_variants( tcr_pose );
	core::scoring::ScoreFunctionOP lowres_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
	core::util::switch_to_residue_type_set( tcr_pose, core::chemical::CENTROID_t );
	//protocols::loops::Loops::const_iterator loop, end;
	protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC ploopmover(cdr3loops, lowres_scorefxn_ );
	ploopmover.apply( tcr_pose );
	core::util::switch_to_residue_type_set( tcr_pose, core::chemical::FULL_ATOM_t );
	//recover_sidechains.apply( tcr_pose );
	repack_tcr_model(tcr_pose);
	protocols::loops::remove_cutpoint_variants( tcr_pose, true );
	tcr_pose.fold_tree( original_ftree );
	//refine cdr3 loops or do post minimization after perturb
	if ( refine ) {
		protocols::loops::loop_mover::refine::LoopMover_Refine_KIC looprefinemover (cdr3loops, core::scoring::get_score_function());
		looprefinemover.apply( tcr_pose );
	} else {
		minimize_tcr_model(tcr_pose);
	}

	return;
}

void ind_remodel_tcr_cdr_loops( core::pose::Pose &tcr_pose, core::Size const loopstart, core::Size const loopend ) {
	//Size num_cycles = option[ max_cycles ];
	core::Size num_cycles = 20;
	core::pose::Pose start_pose = tcr_pose;
	//tcr_pose.dump_pdb("tcr_pose.pdb");
	protocols::loops::Loop cdrloop( loopstart, loopend, 0, 0, false );
	//protocols::loops::remove_cutpoint_variants( tcr_pose, true );
	//protocols::antibody::simple_one_loop_fold_tree( tcr_pose, cdrloop );
	//protocols::loops::add_single_cutpoint_variant( tcr_pose, cdrloop );
	protocols::simple_moves::SwitchResidueTypeSetMover cent( core::chemical::CENTROID );
	cent.apply( tcr_pose );
	protocols::loops::Loop input_loop = cdrloop;
	protocols::antibody::simple_one_loop_fold_tree( tcr_pose, input_loop );
	protocols::loops::LoopsOP pass_loops( new protocols::loops::Loops() );
	pass_loops->add_loop( input_loop );
	pass_loops->set_extended(true);
	core::Real cen_cst_ = 10.0;
	std::string remodel_mover_type = "perturb_kic";
	core::scoring::ScoreFunctionOP lowres_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
	lowres_scorefxn_->set_weight( core::scoring::chainbreak, 10./3. );
	lowres_scorefxn_->set_weight( core::scoring::atom_pair_constraint, cen_cst_ );
	protocols::loops::loop_mover::IndependentLoopMoverOP remodel_mover;
	remodel_mover = utility::pointer::static_pointer_cast< protocols::loops::loop_mover::IndependentLoopMover >
		( protocols::loops::LoopMoverFactory::get_instance()->create_loop_mover(remodel_mover_type, pass_loops) ) ;
	remodel_mover->set_scorefxn( lowres_scorefxn_ );
	core::Size cycle ( 1 );
	while ( cycle < num_cycles ) {
		remodel_mover->apply(tcr_pose);
		++cycle;
	}
	protocols::simple_moves::SwitchResidueTypeSetMover fa( core::chemical::FA_STANDARD );
	fa.apply( tcr_pose );
	utility::vector1<bool> allow_chi_copy( tcr_pose.total_residue(), true );
	for ( core::Size i=loopstart; i<=loopend; i++ ) {
		allow_chi_copy[i] = false;
	}
	protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose, allow_chi_copy );
	recover_sidechains.apply( tcr_pose );
	minimize_tcr_model(tcr_pose);
	return;
}

} //protocols
} //tcr


