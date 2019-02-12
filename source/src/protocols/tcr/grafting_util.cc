// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/grafting_util.cc
/// @brief grafting utility functions for protocols/tcr
/// @author Ragul (ragul@umd.edu)

#include <protocols/grafting/util.hh>
#include <basic/Tracer.hh>
#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/util.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/antibody/util.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <core/pose/chains_util.hh>

static basic::Tracer TR( "protocols.tcr.grafting_util" );

namespace protocols {
namespace tcr {

void graft_cdr_to_fw(core::pose::Pose &scafold, core::pose::Pose const &cdr_piece, core::Size const &start, core::Size const &end, core::Size &nter_overhang, core::Size &cter_overhang ) {
	protocols::grafting::CCDEndsGraftMover ccd_graft_mover;
	ccd_graft_mover = protocols::grafting::CCDEndsGraftMover(start, end, cdr_piece, cter_overhang, nter_overhang);
	ccd_graft_mover.set_cycles(50);
	ccd_graft_mover.final_repack(true);
	ccd_graft_mover.stop_at_closure(true);
	ccd_graft_mover.idealize_insert(true);
	ccd_graft_mover.copy_pdbinfo(true);
	ccd_graft_mover.apply(scafold);
	//Check if graft closed
	//use AnchoredGraftMover if not closed
	using namespace protocols::loops;
	LoopsOP currloop_set( new Loops() );
	Loop Nter_loop = Loop(start-1,start+1,start);
	Loop Cter_loop = Loop(end-1,end+1,end-1);
	currloop_set->add_loop(Nter_loop);
	currloop_set->add_loop(Cter_loop);
	if ( !grafting::graft_closed(scafold, *currloop_set) ) {
		TR<<"CCDEndsGraftMover: Graft not closed"<<std::endl;
		TR<<"Using AnchoredGraftMover"<<std::endl;
		protocols::grafting::AnchoredGraftMover anchored_graft_mover;
		anchored_graft_mover = protocols::grafting::AnchoredGraftMover(start, end, cdr_piece, cter_overhang, nter_overhang);
		anchored_graft_mover.set_cycles(50);
		anchored_graft_mover.final_repack(true);
		anchored_graft_mover.stop_at_closure(true);
		anchored_graft_mover.idealize_insert(true);
		anchored_graft_mover.copy_pdbinfo(true);
		anchored_graft_mover.apply(scafold);
	}
	if ( !grafting::graft_closed(scafold, *currloop_set) ) {
		TR<<"AnchoredGraftMover: Graft not closed"<<std::endl;
	}
	return;
}

void graft_framework(TCRmodel::tmpltinfo &currtmplt, std::string const &template_pdb_path, TCRseqInfo::tcrposi const &aho_pos) {
	core::pose::Pose fr_tmplt;
	core::pose::Pose fr_piece;
	core::Size startcap = aho_pos.cap.begin;
	core::Size endcap = aho_pos.cap.end;
	char tmplt_chainid = currtmplt.tid[5];
	core::import_pose::pose_from_file( fr_tmplt, template_pdb_path, core::import_pose::PDB_file);
	core::Size tmplt_chainnum = get_chain_id_from_chain(tmplt_chainid, fr_tmplt);
	core::pose::Pose fw_tmplt( *fr_tmplt.split_by_chain(tmplt_chainnum) );
	core::Size start_anchor = fw_tmplt.pdb_info()->pdb2pose( tmplt_chainid, aho_pos.cdr1.begin );
	core::Size end_anchor = fw_tmplt.pdb_info()->pdb2pose( tmplt_chainid, aho_pos.cdr3.end );
	fr_piece = protocols::grafting::return_region(fw_tmplt, start_anchor-startcap, end_anchor+endcap);
	currtmplt.tpiece =  core::pose::PoseOP( new core::pose::Pose(fr_piece) );
	return;
}

void graft_cdr(TCRmodel::tmpltinfo &currtmplt, std::string const &template_pdb_path, core::Size const &cdr_start, core::Size const &cdr_end) {
	core::pose::Pose cdr_piece;
	core::import_pose::pose_from_file( cdr_piece, template_pdb_path, core::import_pose::PDB_file);
	char cdr_tmplt_chainid = currtmplt.tid[5];
	core::Size start = cdr_piece.pdb_info()->pdb2pose( cdr_tmplt_chainid, cdr_start );
	core::Size end = cdr_piece.pdb_info()->pdb2pose( cdr_tmplt_chainid, cdr_end );
	core::pose::Pose temp_pose = protocols::grafting::return_region(cdr_piece, start, end);
	currtmplt.tpiece = core::pose::PoseOP( new core::pose::Pose(temp_pose) );
	return;
}

} //protocols
} //tcr


