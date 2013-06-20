// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyInfo.cc
/// @brief Class for getting antibody-specific objects and information
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Project Headers

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/CDRClusterEnum.hh>
#include <protocols/antibody/CDRClusterEnumManager.hh>

// Core Headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/ScoreType.hh>    // scoring stuff can probably be removed / moved elsewhere
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>

// Protocol Headers
#include <protocols/antibody/util.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>  //really?
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/vector1.hh>

//Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Basic headers
#include <numeric/NumericTraits.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <math.h>
#include <utility/string_util.hh>

// Boost headers
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

///////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR("antibody.AntibodyInfo");

namespace protocols {
namespace antibody {





AntibodyInfo::AntibodyInfo( pose::Pose const & pose,
                            AntibodyNumberingSchemeEnum const & numbering_scheme,
                            bool const & cdr_pdb_numbered) :
	framework_info_(NULL)

{
	set_default();

	numbering_scheme_ = numbering_scheme;
	cdr_pdb_numbered_ = cdr_pdb_numbered;

	identify_antibody(pose);
	init(pose);
}

AntibodyInfo::~AntibodyInfo() {}

void AntibodyInfo::set_default() {
	is_camelid_ = false;
	InputPose_has_antigen_ = false;
	predicted_H3_base_type_ = Kinked;
	loopsop_having_allcdrs_=NULL;
	numbering_scheme_ = Aroop;
}

void AntibodyInfo::init(pose::Pose const & pose) {

	if(is_camelid_) total_cdr_loops_ = camelid_last_loop;
	else            total_cdr_loops_ = num_cdr_loops;

	antibody_manager_ = new AntibodyEnumManager();
	cdr_cluster_manager_ = new CDRClusterEnumManager();

	setup_numbering_info_for_scheme(numbering_scheme_);

	setup_CDRsInfo(pose) ;

	setup_FrameWorkInfo(pose) ;

	setup_VL_VH_packing_angle( pose );

	predict_H3_base_type( pose );

}

void AntibodyInfo::setup_numbering_info_for_scheme(AntibodyNumberingSchemeEnum const & numbering_scheme) {

	//////////////////////////////////////////////////////////////////////
	///Feb 2013-
	///Refactored by Jared Adolf-Bryfogle to use Enums and hold the information rather then creating the vectors
	// every time this function was called.
	// Old code is left for reference.


	//Setup the vector

	//TODO: Move this for loop setup stuff to a separate function and then
	//      it can be done in one line:
	//      AntibodyNumbering(num_cdr_loops, vector1< core::Size >(AntibodyNumberingEnum_total, 0))
	//      This wasn't done here because it requires knowing the internal
	//      structure of the AntibodyNumbering typedef, which is not stylish.
	//CDRs
	cdr_numbering_.resize(num_cdr_loops);
	for (core::Size i =1; i <= num_cdr_loops; ++i) {
		cdr_numbering_[i].resize(AntibodyNumberingEnum_total);
	}

	//PackingAngleNumbers
	packing_angle_numbering_.resize(PackingAngleEnum_total);
	for (core::Size i=1; i <= PackingAngleEnum_total; ++i) {
		packing_angle_numbering_[i].resize(AntibodyNumberingEnum_total);
	}


	//**********************************************************************************
	//  Aroop Numbering                                                                *
	//    citation:
	//**********************************************************************************
	if(numbering_scheme == Aroop ) {


		//First and last residue number of each cdr
		cdr_numbering_[h1][start] = 26;   cdr_numbering_[h1][stop] = 35;
		cdr_numbering_[h2][start] = 50;   cdr_numbering_[h2][stop] = 65;
		cdr_numbering_[h3][start] = 95;   cdr_numbering_[h3][stop] = 102;


		cdr_numbering_[l1][start] = 24;    cdr_numbering_[l1][stop] = 34;
		cdr_numbering_[l2][start] = 50;    cdr_numbering_[l2][stop] = 56;
		cdr_numbering_[l3][start] =  89;   cdr_numbering_[l3][stop] = 97;

		// VL-VH packing angle residues
		packing_angle_numbering_[VL_sheet_1][start] = 35;
		packing_angle_numbering_[VL_sheet_1][stop] = 38;
		packing_angle_numbering_[VL_sheet_2][start] = 85;
		packing_angle_numbering_[VL_sheet_2][stop] = 88;

		packing_angle_numbering_[VH_sheet_1][start] = 36;
		packing_angle_numbering_[VH_sheet_1][stop] = 39;
		packing_angle_numbering_[VH_sheet_2][start] = 89;
		packing_angle_numbering_[VH_sheet_2][stop] = 92;
	}
	//**********************************************************************************
	//  Chothia Numbering                                                              *
	//    citation:
	//**********************************************************************************
	else if(numbering_scheme == Chothia ) {

		//First and last residue number of each cdr
		cdr_numbering_[l1][start] = 24;     cdr_numbering_[l1][stop] = 34;
		cdr_numbering_[l2][start] = 50;     cdr_numbering_[l2][stop] = 56;
		cdr_numbering_[l3][start] =  89;    cdr_numbering_[l3][stop] = 97;

		cdr_numbering_[h1][start] = 26;    cdr_numbering_[h1][stop] = 32;
		cdr_numbering_[h2][start] = 52;    cdr_numbering_[h2][stop] = 56;
		cdr_numbering_[h3][start] = 95;    cdr_numbering_[h3][stop] = 102;

		// VL-VH packing angle residues
		packing_angle_numbering_[VL_sheet_1][start] = 35;
		packing_angle_numbering_[VL_sheet_1][stop] = 38;
		packing_angle_numbering_[VL_sheet_2][start] = 85;
		packing_angle_numbering_[VL_sheet_2][stop] = 88;

		packing_angle_numbering_[VH_sheet_1][start] = 36;
		packing_angle_numbering_[VH_sheet_1][stop] = 39;
		packing_angle_numbering_[VH_sheet_2][start] = 89;
		packing_angle_numbering_[VH_sheet_2][stop] = 92;

		// VL-VH packing angle residues
		//pack_angle_start.push_back(35); pack_angle_stop.push_back(38);  //VL
		//pack_angle_start.push_back(85); pack_angle_stop.push_back(88);  //VL
		//pack_angle_start.push_back(36); pack_angle_stop.push_back(39);  //VH
		//pack_angle_start.push_back(89); pack_angle_stop.push_back(92);  //VH
	}

	//**********************************************************************************
	//  Kabat Numbering                                                                *
	//    citation:
	//**********************************************************************************
	else if(numbering_scheme == Kabat ) {
	}
	//**********************************************************************************
	//  Enhanced_Chothia Numbering                                                     *
	//    Abhinandan, K.R. and Martin, A.C.R. (2008) Immunology, 45, 3832-3839.        *
	//**********************************************************************************
	else if(numbering_scheme == Enhanced_Chothia) {
	}
	//**********************************************************************************
	//  AHO Numbering                                                                  *
	//    A. Honegger & A. Plückthun. (2001) J. Mol. Biol, 309 (2001)657-670
	//**********************************************************************************
	else if(numbering_scheme == AHO) {
	}
	//**********************************************************************************
	//  Modified Honegger-Plückthun (Aho) Numbering                                                                   *
	//    North, B., A. Lehmann, et al. (2011). JMB 406(2): 228-256. (Not described in paper)
	//**********************************************************************************
	else if(numbering_scheme == Modified_AHO) {

		//First and last residue number of each cdr
		cdr_numbering_[l1][start] = 24;     cdr_numbering_[l1][stop] = 42;
		cdr_numbering_[l2][start] = 57;     cdr_numbering_[l2][stop] = 72;
		cdr_numbering_[l3][start] =  107;  cdr_numbering_[l3][stop] = 138;

		cdr_numbering_[h1][start] = 24;    cdr_numbering_[h1][stop] = 42;
		cdr_numbering_[h2][start] = 57;    cdr_numbering_[h2][stop] = 69;
		cdr_numbering_[h3][start] = 107;  cdr_numbering_[h3][stop] = 138;

		// VL-VH packing angle residues - Look at the symmetry!
		packing_angle_numbering_[VL_sheet_1][start] = 43;
		packing_angle_numbering_[VL_sheet_1][stop] = 46;
		packing_angle_numbering_[VL_sheet_2][start] = 103;
		packing_angle_numbering_[VL_sheet_2][stop] = 106;

		packing_angle_numbering_[VH_sheet_1][start] = 43;
		packing_angle_numbering_[VH_sheet_1][stop] = 46;
		packing_angle_numbering_[VH_sheet_2][start] = 103;
		packing_angle_numbering_[VH_sheet_2][stop] = 106;
	}
	//**********************************************************************************
	//  IMGT Numbering                                                                 *
	//    citation:
	//**********************************************************************************
	else if(numbering_scheme == IMGT) {
	} else {
		utility_exit_with_message("the numbering schemes can only be 'Aroop','Chothia','Kabat', 'Enhanced_Chothia', 'AHO', Modified_AHO', 'IMGT' !!!!!! ");
	}

	//local_numbering_info.push_back(start);
	//local_numbering_info.push_back(stop);
	//local_numbering_info.push_back(pack_angle_start);
	//local_numbering_info.push_back(pack_angle_stop);
	//return local_numbering_info;
}

std::string
AntibodyInfo::get_Current_AntibodyNumberingScheme() {
	return antibody_manager_->numbering_scheme_enum_to_string(numbering_scheme_);
}

void AntibodyInfo::identify_antibody(pose::Pose const & pose) {
	//Time to rewrite.  Who needs to have L and H as the first chains? come on.

	//Jadofbr (4/2013)  Allow any order in the PDB file by adding this for loop..
	//Todo: Add command-line argument for other types (SCFv, nanobody, Diabody, etc.)
	bool H_found = false;
	vector1<char> H_sequence;
	bool L_found = false;
	vector1<char> L_sequence;
	for (core::Size i = 1; i<= pose.conformation().num_chains(); ++i) {
		if(   pose.pdb_info()->chain(pose.conformation().chain_end(i)) == 'L') {
			L_found=true;
			L_chain_ = i;

			for (core::Size x = 1; x<=pose.total_residue(); ++x) {
				//No good way to get a chains sequence that I know about.
				if (pose.residue(x).chain()==L_chain_) {
					L_sequence.push_back(pose.residue(i).name1());
					sequence_map_[x]= pose.residue(x).name1();
				}
			}
		} else if(pose.pdb_info()->chain(pose.conformation().chain_end(i)) == 'H' ) {
			H_found=true;
			H_chain_=i;

			for (core::Size x = 1; x<=pose.total_residue(); ++x) {
				//No good way to get a chains sequence that I know about.
				if (pose.residue(x).chain()==H_chain_) {
					H_sequence.push_back(pose.residue(i).name1());
					sequence_map_[x]= pose.residue(x).name1();
				}
			}
		} else {
			continue;
		}
	}


	switch (pose.conformation().num_chains() ) {
	case 0:
		throw excn::EXCN_Msg_Exception("the number of chains in the input pose is '0' !!");
		break;
	case 1:
		//if pose has only "1" chain, it is a nanobody
		if ( H_found ) {
			is_camelid_ = true;
			InputPose_has_antigen_=false;
		} else {
			throw excn::EXCN_Msg_Exception("  A): the input pose has only 1 chain, if it is a nanobody, the chain ID is supposed to be 'H' !!");
		}
		break;
	case 2:
		// if pose has "2" chains, it can be 2 possibilities
		// possiblity 1): L and H, regular antibody
		if (  L_found  && H_found ) {
			is_camelid_ = false;
			InputPose_has_antigen_ = false;
		}
		// possiblity 2): H nanobody and antigen
		else if ( H_found) {
			is_camelid_ = true;
			InputPose_has_antigen_ = true;
		} else {
			throw excn::EXCN_Msg_Exception("  B): the input pose has two chains, 1). if it is nanobody, the chain should be 'H'. 2). Light chain SCFv not implemented ");
		}
		break;
	default:
		// if pose has >=3 chains, it can be 2 possibilities
		// possiblity 1): L and H, and antigen
		if(   L_found  &&  H_found ) {
			is_camelid_ = false;
			InputPose_has_antigen_ = true;
		}
		// possiblity 2):  H nanobody and antigen
		else if ( H_found ) {
			is_camelid_ = true;
			InputPose_has_antigen_ = true;
		} else {
			throw excn::EXCN_Msg_Exception("  C). the input pose has more than two chains, but either 1) Light chain SCFv not implemented 2) Antibody is not renumbered ");
		}
		break;
	}

	/// record the antibody sequence
	if (is_camelid_) {
		ab_sequence_ = H_sequence;
	} else {
		//ab_sequence_.reserve(L_sequence.size() + H_sequence.size());
		ab_sequence_.insert(ab_sequence_.end(), L_sequence.begin(), L_sequence.end());
		ab_sequence_.insert(ab_sequence_.end(), H_sequence.begin(), H_sequence.end());

	}
}



//Jadolfbr 3/2012  Refactored to be numbering scheme agnostic.
void AntibodyInfo::setup_CDRsInfo( pose::Pose const & pose ) {


	for (Size i=1; i<=3; ++i) {
		Chain_IDs_for_CDRs_.push_back('H');    // HEAVY chain first
	}
	for (Size i=1; i<=3; ++i) {
		Chain_IDs_for_CDRs_.push_back('L');    // light
	}

	int loop_start_in_pose, loop_stop_in_pose, cut_position ;
	loopsop_having_allcdrs_ = new loops::Loops();

	for (Size i=start_cdr_loop; i<=total_cdr_loops_; ++i ) {
		loop_start_in_pose = pose.pdb_info()->pdb2pose( Chain_IDs_for_CDRs_[i], cdr_numbering_[i][start]);
		if(i != h3 ) {
			loop_stop_in_pose= pose.pdb_info()->pdb2pose( Chain_IDs_for_CDRs_[i], cdr_numbering_[i][stop]);
			cut_position = (loop_stop_in_pose - loop_start_in_pose +1) /2 + loop_start_in_pose;
		} else {
			loop_stop_in_pose  = pose.pdb_info()->pdb2pose( Chain_IDs_for_CDRs_[i], cdr_numbering_[i][stop]+1 );
			loop_stop_in_pose -=1;
			// JQX:
			// One should always see 95-102 as the positions for your H3 in your FR02.pdb, but as a matter of fact,
			// the antibody script just copied h3.pdb (heavy atoms) into the FR02.pdb, sometimes one sees the stop
			// postition pdb number 98, not 102, if the h3.pdb is short. Therefore, one useing the pdb number 102 to
			// define h3 fails!
			// of PDB number 103, then minus 1 will give you the last residue of h3.
			cut_position = (loop_start_in_pose +1 ) ;
			// JQX:
			// why this is different compared to other cuts of other loops?
			// Aroop seems did this in his old R3 code, CHECK LATER !!!
		}

		loops::Loop  one_loop(loop_start_in_pose, loop_stop_in_pose, cut_position);
		loops::LoopsOP one_loops = new loops::Loops();
		one_loops->add_loop(one_loop);

		// make a "LoopsOP" object, in which each "Loop" was saved
		loopsop_having_allcdrs_->add_loop(one_loop);

		// make a "vector1" of "LoopsOP" object, each "LoopsOP" has one "Loop" object
		vector1_loopsop_having_cdr_.push_back(one_loops);
	}

	/// FIXME:  ***********************
	loopsop_having_allcdrs_->sequential_order(); /// TODO: kind of dangerous here

	TR<<"Successfully finished the CDR definition"<<std::endl;

}

vector1< vector1<FrameWork> >
AntibodyInfo::get_AntibodyFrameworkInfo() const {
	if (framework_info_.empty()) {
		utility_exit_with_message("Numbering scheme setup failed for Framework.");
	} else if (numbering_scheme_ != Aroop) {
		utility_exit_with_message("Framework info only currently works for Aroop numbering scheme.");
	} else {
		return framework_info_;
	}
}

H3BaseTypeEnum
AntibodyInfo::get_Predicted_H3BaseType() const {
	if (numbering_scheme_ != Aroop) {
		utility_exit_with_message("H3 Base type is only currently correct for Aroop numbering scheme.");
	} else {
		return predicted_H3_base_type_;
	}
}

void AntibodyInfo::setup_FrameWorkInfo( pose::Pose const & pose ) {

	//TODO jadolfbr.  This is a hardcoded nightmare.
	//Used only by align_to_native in antibody_util.  Needs to be refactored for other numbering schemes

	FrameWork frmwk;
	vector1<FrameWork> Lfr, Hfr;

	core::Size H_begin_pos_num = 0;
	core::Size H_end_pos_num = 0;
	core::Size L_begin_pos_num = 0;
	core::Size L_end_pos_num = 0;


	if(is_camelid() == true ) {
		H_begin_pos_num=pose.conformation().chain_begin(H_chain_);
		H_end_pos_num=pose.conformation().chain_end(H_chain_);
	} else {
		L_begin_pos_num=pose.conformation().chain_begin(L_chain_);
		L_end_pos_num=pose.conformation().chain_end(L_chain_);
		H_begin_pos_num=pose.conformation().chain_begin(H_chain_);
		H_end_pos_num=pose.conformation().chain_end(H_chain_);;


		if (  L_begin_pos_num   >=    pose.pdb_info()->pdb2pose('L',23)   )  {
			throw excn::EXCN_Msg_Exception( "L chain 1st residue starting after L 23, framework definition failed!!! " );
		}
		if (  L_end_pos_num     <=    pose.pdb_info()->pdb2pose('L', 97)    ) {
			throw excn::EXCN_Msg_Exception( "L chain last residue ending before L 97, framework definition failed!!! " );
		}
	}


	if (    H_begin_pos_num    >=     pose.pdb_info()->pdb2pose('H', 26 )      ) {
		throw excn::EXCN_Msg_Exception( "H chain 1st residue starting after H 26, framework definition failed!!! " );
	}

	if (    H_end_pos_num      <=     pose.pdb_info()->pdb2pose('H', 103)      ) {
		throw excn::EXCN_Msg_Exception( "H chain last residue ending before H 103, framework definition failed!!! " );
	}


	switch (numbering_scheme_) {
	case Aroop:
		if(! is_camelid_) {
			frmwk.chain_name='L';
			if(   L_begin_pos_num <= pose.pdb_info()->pdb2pose('L',5)      ) { // <= 5
				frmwk.start=pose.pdb_info()->pdb2pose('L',5);
				frmwk.stop=pose.pdb_info()->pdb2pose('L',6);
				Lfr.push_back(frmwk);
			} else if (     L_begin_pos_num   <=    pose.pdb_info()->pdb2pose('L',6)        ) { // 5 <= x <= 6
				frmwk.start=L_begin_pos_num;
				frmwk.stop=pose.pdb_info()->pdb2pose('L',6);
				Lfr.push_back(frmwk);
			}
			if(   L_begin_pos_num   <=    pose.pdb_info()->pdb2pose('L',10)      ) { // <= 10
				frmwk.start=pose.pdb_info()->pdb2pose('L',10);
				frmwk.stop=pose.pdb_info()->pdb2pose('L',23);
				Lfr.push_back(frmwk);
			} else if (   L_begin_pos_num   <=    pose.pdb_info()->pdb2pose('L',23)           ) { //  10 <= x <=23
				frmwk.start=L_begin_pos_num;
				frmwk.stop=pose.pdb_info()->pdb2pose('L',23);
				Lfr.push_back(frmwk);
			}
			frmwk.start=pose.pdb_info()->pdb2pose('L',35);
			frmwk.stop=pose.pdb_info()->pdb2pose('L',38);
			Lfr.push_back(frmwk);
			frmwk.start=pose.pdb_info()->pdb2pose('L',45);
			frmwk.stop=pose.pdb_info()->pdb2pose('L',49);
			Lfr.push_back(frmwk);
			frmwk.start=pose.pdb_info()->pdb2pose('L',57);
			frmwk.stop=pose.pdb_info()->pdb2pose('L',66);
			Lfr.push_back(frmwk);
			frmwk.start=pose.pdb_info()->pdb2pose('L',71);
			frmwk.stop=pose.pdb_info()->pdb2pose('L',88);
			Lfr.push_back(frmwk);
			if (   L_end_pos_num   >=    pose.pdb_info()->pdb2pose('L',105)      ) {
				frmwk.start=pose.pdb_info()->pdb2pose('L',98);
				frmwk.stop=pose.pdb_info()->pdb2pose('L',105);
				Lfr.push_back(frmwk);
			} else {
				frmwk.start=pose.pdb_info()->pdb2pose('L',98);
				frmwk.stop=L_end_pos_num;
				Lfr.push_back(frmwk);
			}
		}

		frmwk.chain_name='H';
		if(   H_begin_pos_num   <=    pose.pdb_info()->pdb2pose('H',5)      ) { // <= 5
			frmwk.start=pose.pdb_info()->pdb2pose('H',5);
			frmwk.stop=pose.pdb_info()->pdb2pose('H',6);
			Hfr.push_back(frmwk);
		} else if (   H_begin_pos_num   <=    pose.pdb_info()->pdb2pose('H',6)      ) { // 5 <= x <= 6
			frmwk.start=H_begin_pos_num;
			frmwk.stop=pose.pdb_info()->pdb2pose('H',6);
			Hfr.push_back(frmwk);
		}
		if(   H_begin_pos_num   <=    pose.pdb_info()->pdb2pose('H',10)      ) { // <= 10
			frmwk.start=pose.pdb_info()->pdb2pose('H',10);
			frmwk.stop=pose.pdb_info()->pdb2pose('H',25);
			Hfr.push_back(frmwk);
		} else if(   H_begin_pos_num   <=    pose.pdb_info()->pdb2pose('H',25)      ) { //  10 <= x <=25
			frmwk.start=H_begin_pos_num;
			frmwk.stop=pose.pdb_info()->pdb2pose('H',25);
			Hfr.push_back(frmwk);
		}
		frmwk.start=pose.pdb_info()->pdb2pose('H',36);
		frmwk.stop=pose.pdb_info()->pdb2pose('H',39);
		Hfr.push_back(frmwk);
		frmwk.start=pose.pdb_info()->pdb2pose('H',46);
		frmwk.stop=pose.pdb_info()->pdb2pose('H',49);
		Hfr.push_back(frmwk);
		frmwk.start=pose.pdb_info()->pdb2pose('H',66);
		frmwk.stop=pose.pdb_info()->pdb2pose('H',94);
		Hfr.push_back(frmwk);
		if(   H_end_pos_num >=  pose.pdb_info()->pdb2pose('H',110)         ) {
			frmwk.start=pose.pdb_info()->pdb2pose('H',103);
			frmwk.stop=pose.pdb_info()->pdb2pose('H',110);
			Hfr.push_back(frmwk);
		} else {
			frmwk.start=pose.pdb_info()->pdb2pose('H',103);
			frmwk.stop=H_end_pos_num;
			Hfr.push_back(frmwk);
		}
		break;
	case Chothia:
		break;
	case Kabat:
		break;
	case Enhanced_Chothia:
		break;
	case AHO:
		break;
	case Modified_AHO:
		break;
	case IMGT:
		break;
	default:
		utility_exit_with_message("the numbering schemes can only be 'Aroop','Chothia','Kabat', 'Enhanced_Chothia', 'AHO', 'IMGT' !!!!!! ");
		break;
	}

	if (Lfr.size()>0)  {
		framework_info_.push_back(Lfr);
	}
	if (Hfr.size()>0)  {
		framework_info_.push_back(Hfr);
	}

	//Removed.  Check is now in get_AntibodyFrameworkInfo.
	//else {utility_exit_with_message("The heavy chain has no framework? This cannot be correct");}


}


void AntibodyInfo::setup_VL_VH_packing_angle( pose::Pose const & pose ) {

	vector1<char> Chain_IDs_for_packing_angle;
	Chain_IDs_for_packing_angle.resize(PackingAngleEnum_total);
	Chain_IDs_for_packing_angle[VL_sheet_1] = 'L';
	Chain_IDs_for_packing_angle[VL_sheet_2] = 'L';
	Chain_IDs_for_packing_angle[VH_sheet_1] = 'H';
	Chain_IDs_for_packing_angle[VH_sheet_2] = 'H';

	Size packing_angle_start_in_pose, packing_angle_stop_in_pose;

	for (Size i=1; i<=PackingAngleEnum_total; ++i) {
		packing_angle_start_in_pose = pose.pdb_info()->pdb2pose( Chain_IDs_for_packing_angle[i], packing_angle_numbering_[i][start]);
		packing_angle_stop_in_pose = pose.pdb_info()->pdb2pose( Chain_IDs_for_packing_angle[i], packing_angle_numbering_[i][stop]);
		for (Size j=packing_angle_start_in_pose; j<=packing_angle_stop_in_pose; j++) {
			packing_angle_residues_.push_back( j );
		}
	}
}




////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
///	predicit H3 cterminus base type (Kinked or Extended) based on sequence   ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
void AntibodyInfo::predict_H3_base_type( pose::Pose const & pose ) {
	if( is_camelid_ ) {
		detect_and_set_camelid_CDR_H3_stem_type( pose );
	} else {
		//detect_and_set_regular_CDR_H3_stem_type( pose );
		detect_and_set_regular_CDR_H3_stem_type_new_rule( pose );
	}
} // detect_CDR_H3_stem_type


void AntibodyInfo::detect_and_set_camelid_CDR_H3_stem_type(pose::Pose const & pose ) {
	TR << "AC Detecting Camelid CDR H3 Stem Type" << std::endl;

	bool kinked_H3 (false);
	bool extended_H3 (false);

	// extract single letter aa codes for the chopped loop residues
	vector1< char > cdr_h3_sequence;
	for( Size ii = get_CDR_loop(h3).start() - 2; ii <= get_CDR_loop(h3).stop(); ++ii )
		cdr_h3_sequence.push_back( pose.sequence()[ii-1] );

	// Rule for extended
	if( ( ( get_CDR_loop(h3).stop() - get_CDR_loop(h3).start() ) ) >= 12 ) {
		if( ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'Y' ) ||
		        ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'W' ) ||
		        ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'F' ) ) &&
		        ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'H' ) &&
		        ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'G' ) )
			extended_H3 = true;
	}

	if( !extended_H3 ) {
		kinked_H3 = true;
		if(           ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'R' ) ||
		              ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'Y' ) ||
		              ((     ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'W' )    ) &&
		               (     ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'W' )    ) &&
		               ( 	  ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'W' )    ))
		  )
			kinked_H3 = false;
	}

	if (kinked_H3) predicted_H3_base_type_ = Kinked;
	if (extended_H3) predicted_H3_base_type_ = Extended;
	if (!kinked_H3 && !extended_H3) predicted_H3_base_type_ = Neutral;
	TR << "AC Finished Detecting Camelid CDR H3 Stem Type: " << antibody_manager_->h3_base_type_enum_to_string(predicted_H3_base_type_) << std::endl;
}


void AntibodyInfo::detect_and_set_regular_CDR_H3_stem_type( pose::Pose const & pose ) {
	TR << "AC Detecting Regular CDR H3 Stem Type" << std::endl;
	bool extended_H3 (false) ;
	bool kinked_H3 (false);
	bool is_H3( false );

	// extract single letter aa codes for the chopped loop residues
	vector1< char > cdr_h3_sequence;
	for( Size ii = get_CDR_loop(h3).start() - 2; ii <= get_CDR_loop(h3).stop()+1 ; ++ii )
		cdr_h3_sequence.push_back( pose.sequence()[ii-1] );

	// Rule 1a for standard kink
	if( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'D') {
		kinked_H3 = true;
		is_H3 = true;
	}

	// Rule 1b for standard extended form
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D')
	        && ( (cdr_h3_sequence[2] != 'K') &&
	             (cdr_h3_sequence[2] != 'R') ) && (is_H3 != true)) {
		extended_H3 = true;
		is_H3 = true;
	}

	if( !is_H3 ) {
		// Rule 1b extension for special kinked form
		bool is_basic( false ); // Special basic residue exception flag
		for(Size ii = 3; ii <= Size(cdr_h3_sequence.size() - 4); ++ii) {
			if( cdr_h3_sequence[ii] == 'R' || cdr_h3_sequence[ii] == 'K') {
				is_basic = true;
				break;
			}
		}

		if( !is_basic ) {
			Size L49_pose_number = pose.pdb_info()->pdb2pose( 'L', 49 );
			char aa_code_L49 = pose.residue( L49_pose_number ).name1();
			if( aa_code_L49 == 'R' || aa_code_L49 == 'K')
				is_basic = true;
		}
		if( is_basic ) {
			kinked_H3 = true;
			is_H3 = true;
		}
	}

	// Rule 1c for kinked form with salt bridge
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
	        ( (cdr_h3_sequence[2] == 'K') ||
	          (cdr_h3_sequence[2] == 'R') ) &&
	        ( (cdr_h3_sequence[1] != 'K') &&
	          (cdr_h3_sequence[1] != 'R') ) && (is_H3 != true) ) {
		kinked_H3 = true;
		is_H3 = true;
		if( !is_H3 ) {
			bool is_basic( false ); // Special basic residue exception flag
			Size L46_pose_number = pose.pdb_info()->pdb2pose( 'L', 46 );
			char aa_code_L46 = pose.residue( L46_pose_number ).name1();
			if( aa_code_L46 == 'R' || aa_code_L46 == 'K')
				is_basic = true;
			if( is_basic ) {
				extended_H3 = true;
				is_H3 = true;
			}
		}
	}

	// Rule 1d for extened form with salt bridge
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
	        ( ( cdr_h3_sequence[ 2 ] == 'K') ||
	          (cdr_h3_sequence[2] == 'R')) &&
	        ( (cdr_h3_sequence[1] == 'K') ||
	          (cdr_h3_sequence[1] == 'R') ) && (is_H3 != true) ) {
		extended_H3 = true;
		is_H3 = true;
	}

	if (kinked_H3) predicted_H3_base_type_ = Kinked;
	if (extended_H3) predicted_H3_base_type_ = Extended;
	if (!kinked_H3 && !extended_H3) predicted_H3_base_type_ = Neutral;
	TR << "AC Finished Detecting Regular CDR H3 Stem Type: " << antibody_manager_->h3_base_type_enum_to_string(predicted_H3_base_type_) << std::endl;
} // detect_regular_CDR_H3_stem_type()


void AntibodyInfo::detect_and_set_regular_CDR_H3_stem_type_new_rule( pose::Pose const & pose ) {
	TR << "AC Detecting Regular CDR H3 Stem Type" << std::endl;

	if (numbering_scheme_ != Aroop) {
		TR << "Stem type could not be set."<<std::endl;
		predicted_H3_base_type_ = Unknown;
		return;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////HACK/////////////////////////////////////////////////////


	bool extended_H3 (false) ;
	bool kinked_H3 (false);

	// extract single letter aa codes for the chopped loop residues
	vector1< char > cdr_h3_sequence;
	for( Size ii = get_CDR_loop(h3).start() - 2; ii <= get_CDR_loop(h3).stop() + 1; ++ii )
		cdr_h3_sequence.push_back( pose.sequence()[ii-1] );
	//for (Size i=1; i<=cdr_h3_sequence.size();++i){    TR<<cdr_h3_sequence[i];} TR<<std::endl;

	/// @author: Daisuke Kuroda (dkuroda1981@gmail.com) 06/18/2012
	///
	/// @last_modified 06/18/2012
	///
	/// @reference Kuroda et al. Proteins. 2008 Nov 15;73(3):608-20.
	///			   Koliansnikov et al. J Bioinform Comput Biol. 2006 Apr;4(2):415-24.

	// This is only for rule 1b
	bool is_basic( false ); // Special basic residue exception flag
	if( !is_basic ) {
		Size L49_pose_number = pose.pdb_info()->pdb2pose( 'L', 49 );
		char aa_code_L49 = pose.residue( L49_pose_number ).name1();
		if( aa_code_L49 == 'R' || aa_code_L49 == 'K')
			is_basic = true;
	}

	/// START H3-RULE 2007
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
	        ( ( cdr_h3_sequence[ 2 ] == 'K') || (cdr_h3_sequence[2] == 'R') ) &&
	        ( ( cdr_h3_sequence[ 1 ] == 'K') || (cdr_h3_sequence[1] == 'R') ) ) {
		// Rule 1d for extened form with salt bridge
		extended_H3 = true;
	} else if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
	           ( ( cdr_h3_sequence[ 2 ] == 'K') || ( cdr_h3_sequence[ 2 ] == 'R') ) &&
	           ( ( cdr_h3_sequence[ 1 ] != 'K') && ( cdr_h3_sequence[ 1 ] != 'R') ) ) {
		// Rule 1c for kinked form with salt bridge with/without Notable signal (L46)
		// Special basic residue exception flag
		Size L46_pose_number = pose.pdb_info()->pdb2pose( 'L', 46 );
		char aa_code_L46 = pose.residue( L46_pose_number ).name1();

		// Special Tyr residue exception flag
		Size L36_pose_number = pose.pdb_info()->pdb2pose( 'L', 36 );
		char aa_code_L36 = pose.residue( L36_pose_number ).name1();

		if( ( aa_code_L46 == 'R' || aa_code_L46 == 'K') && aa_code_L36 != 'Y' ) {
			extended_H3 = true;
		} else {
			kinked_H3   = true;
		}
	} else if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D' ) &&
	           ( cdr_h3_sequence[ 2 ] != 'K' ) && ( cdr_h3_sequence[ 2 ] != 'R' ) &&
	           ( is_basic == true ) ) {
		// Rule 1b for standard extended form with Notable signal (L49)
		kinked_H3   = true;
	} else if( ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'F' ) &&
	             ( cdr_h3_sequence[ cdr_h3_sequence.size() - 4 ] == 'A' ) ) ||
	           ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'F' ) &&
	             ( cdr_h3_sequence[ cdr_h3_sequence.size() - 4 ] == 'G' ) ) ||
	           ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'M' ) &&
	             ( cdr_h3_sequence[ cdr_h3_sequence.size() - 4 ] == 'A' ) ) ||
	           ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'M' ) &&
	             ( cdr_h3_sequence[ cdr_h3_sequence.size() - 4 ] == 'G' ) ) ) {
		// This is new feature
		kinked_H3   = true;
	} else if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'R' ) ||
	           ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'K' ) ||
	           ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'D' ) ||
	           ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'N' ) ) {
		// This is new feature
		extended_H3 = true;
	} else if( ( ( cdr_h3_sequence[ 3 ] == 'Y' ) &&
	             ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'F' ) ) ||
	           ( (cdr_h3_sequence[ 3 ] == 'Y' ) &&
	             (cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'M') ) ) {
		// This is new feature
		extended_H3 = true;
	} else if( cdr_h3_sequence.size() - 3  == 7 ) {
		// This is new feature
		extended_H3 = true;
	} else if( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D' ) {
		// Rule 1b for standard extended form without Notable signal (L49)
		extended_H3 = true;
	} else if( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'D' ) {
		// Rule 1a for standard kink. i.e. No sequence feature...
		kinked_H3 = true;
	}
	// END H3-RULE 2007

	if (kinked_H3) predicted_H3_base_type_ = Kinked;
	if (extended_H3) predicted_H3_base_type_ = Extended;
	if (!kinked_H3 && !extended_H3) predicted_H3_base_type_ = Neutral;
	TR << "AC Finished Detecting Regular CDR H3 Stem Type: " << antibody_manager_->h3_base_type_enum_to_string(predicted_H3_base_type_) << std::endl;

	TR << "AC Finished Detecting Regular CDR H3 Stem Type: "
	   << "Kink: " << kinked_H3 << " Extended: " << extended_H3 << std::endl;
} // detect_regular_CDR_H3_stem_type()

std::pair <CDRClusterEnum, Real>
AntibodyInfo::get_CDR_cluster(CDRNameEnum const cdr_name) {
	if (cdr_clusters_.empty()) {
		utility_exit_with_message("CDR clusters not setup.  Cannot proceed.");
	}
	return std::make_pair(cdr_clusters_[cdr_name], cdr_cluster_distances_[cdr_name]);
}

Size
AntibodyInfo::get_CDR_start(CDRNameEnum const & cdr_name, pose::Pose const & pose) const {
	return pose.pdb_info()->pdb2pose(Chain_IDs_for_CDRs_[cdr_name], cdr_numbering_[cdr_name][start]);
}

Size
AntibodyInfo::get_CDR_end(CDRNameEnum const & cdr_name, pose::Pose const & pose) const {
	return pose.pdb_info()->pdb2pose(Chain_IDs_for_CDRs_[cdr_name], cdr_numbering_[cdr_name][stop]);
}


std::vector<Vector>
AntibodyInfo::kink_anion_atoms(const core::pose::Pose & pose) const {
	Size resi = get_CDR_loop(h3).stop() - 1;
	core::conformation::Residue res = pose.residue(resi);
	//print "H3_N-1 (%i): %s" % (resi,res.name3())
	std::vector<Vector> atoms;
	switch (res.name1()) {
	case 'D' :
		atoms.push_back(res.xyz("OD1"));
		atoms.push_back(res.xyz("OD2"));
		break;
	case 'E' :
		atoms.push_back(res.xyz("OD1"));
		atoms.push_back(res.xyz("OD2"));
		break;
	}
	return atoms;
}

/// @brief return side chain cation atoms (typically Lys/His nitrogens) in the kink bulge HBond
std::vector<Vector>
AntibodyInfo::kink_cation_atoms(const core::pose::Pose & pose) const {
	Size resi = get_CDR_loop(h3).start() - 1;
	core::conformation::Residue res = pose.residue(resi);
	//print "H3_0   (%i): %s" % (resi,res.name3())
	std::vector<Vector> atoms;
	switch (res.name1()) {
	case 'R' :
		atoms.push_back(res.xyz("NH1"));
		atoms.push_back(res.xyz("NH2"));
		break;
	case 'K' :
		atoms.push_back(res.xyz("NZ"));
		break;
	}
	return atoms;
}



void
AntibodyInfo::setup_CDR_clusters(pose::Pose const & pose) {

	TR << "Setting up CDR Clusters" << std::endl;
	cdr_clusters_.resize(CDRNameEnum_total);
	cdr_cluster_distances_.resize(CDRNameEnum_total);
	for (core::Size i=1; i<=CDRNameEnum_total; ++i) {
		CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
		set_cdr_cluster(pose, cdr_name);
	}
}

bool
AntibodyInfo::clusters_setup() {
	if (!cdr_clusters_.empty()) {
		return true;
	} else {
		return false;
	}
}

void
AntibodyInfo::set_cdr_cluster(pose::Pose const & pose, CDRNameEnum const & cdr_name) {
	using std::string;
	using std::map;
	using core::Size;
	using core::Real;

	string cdr;
	Size length;
	string cluster;
	string type;
	string fullcluster;
	string ss;
	string phis;
	string psis;
	string omegas;
	Real PI = numeric::NumericTraits<Real>::pi();
	string const path = "sampling/antibodies/cluster_center_dihedrals.txt";


	Size cdr_length = get_CDR_length(cdr_name);
	loops::Loop loop = get_CDR_loop(cdr_name);
	Size start = loop.start();

	//Get current cdr phi, psi, omega Once.
	utility::io::izstream clus_stream;
	basic::database::open(clus_stream, path);
	map <string, vector1< Real > > pose_angles;
	string model_ss;
	Real cis_cutoff = 90.00;
	for (Size resnum = start; resnum<=(start+cdr_length-1); resnum++) {
		pose_angles["phi"].push_back(pose.phi(resnum));
		pose_angles["psi"].push_back(pose.psi(resnum));
		//Pythonic way of doing it, may need to change.
		if (std::abs(pose.omega(resnum))>=cis_cutoff) {
			model_ss = model_ss+"T";
		} else {
			model_ss = model_ss+"C";
		}
	}

	TR <<"CDR: "<<get_CDR_Name(cdr_name)<<" "<<cdr_length<<" "<<model_ss<<std::endl;

	map <Real, string > k_distances_to_cluster;
	vector1< Real > k_distances;
	Size line_count = 0;
	bool cluster_found = false;
	while ( ! clus_stream.eof() ) {
		//TR << "Reading file..line"<<line_count<<std::endl;
		++line_count;
		clus_stream >> cdr >> length >> cluster >> type >> fullcluster >> ss >> phis >> psis >> omegas;
		//Should I do the measuring here, or actually put the data in ram, then go through the data?
		//TR << cdr <<" "<<length<<" "<<ss<<std::endl;
		if (cdr==get_CDR_Name(cdr_name) && length==cdr_length && ss==model_ss) {
			cluster_found=true;
			Real k_distance_to_cluster=0.0;
			//Calculate, add into cluster_distances.
			map <string, vector1 <string> > cluster_angles;

			boost::split(cluster_angles["phi"], phis, boost::is_any_of(","));
			boost::split(cluster_angles["psi"], psis, boost::is_any_of(","));
			//boost::split(cluster_angles["omega"], omegas, boost::is_any_of(":"));

			//Calculate:
			for (Size i=1; i<=cdr_length; i++) {
				//Need to convert angles to positive values for now.
				std::istringstream phi_stream(cluster_angles["phi"][i]);
				std::istringstream psi_stream(cluster_angles["psi"][i]);
				Real phi;
				Real psi;
				phi_stream >> phi;
				psi_stream >> psi;

				Real phi_d = (2 * (1- cos ((pose_angles["phi"][i]-phi)*PI/180)));
				Real psi_d = (2 * (1- cos ((pose_angles["psi"][i]-psi)*PI/180)));
				k_distance_to_cluster = k_distance_to_cluster+phi_d+psi_d;

			}
			k_distances_to_cluster[k_distance_to_cluster]=fullcluster;
			k_distances.push_back(k_distance_to_cluster);
		} else {
			continue;
		}

	}//End while Loop

	clus_stream.close();

	///Take the minimum distance as the cluster.
	std::pair<CDRClusterEnum, Real> result;

	if (! cluster_found) {
		TR <<"Cluster not found for CDR "<<get_CDR_Name(cdr_name)<<std::endl;
		cdr_clusters_[cdr_name] = NA;
		cdr_cluster_distances_[cdr_name] = 1000;
		return;
	} else {
		//Get minimum and set cluster.
		Real d = k_distances[1];
		for (Size i=2; i<k_distances.size(); i++) {
			if (k_distances[i]<d) {
				d = k_distances[i];
			}
		}
		string found_cluster = k_distances_to_cluster[d];
		TR << get_CDR_Name(cdr_name) <<" cluster found as " << found_cluster << " at k_distance: " << d <<std::endl;
		//TR <<"Setting this as closest cluster, as no cutoff is yet set. PLEASE manually compare structures.  " <<std::endl;
		CDRClusterEnum cluster_enum = cdr_cluster_manager_->cdr_cluster_string_to_enum(found_cluster);
		cdr_clusters_[cdr_name] = cluster_enum;
		cdr_cluster_distances_[cdr_name] = d;
		return;
	}
}

std::string
AntibodyInfo::get_CDR_Name(CDRNameEnum const & cdr_name) const {
	return antibody_manager_->cdr_name_enum_to_string(cdr_name);
}

Size
AntibodyInfo::get_CDR_length(CDRNameEnum const & cdr_name) const {
	loops::Loop cdr_loop = get_CDR_loop(cdr_name);
	Size l = cdr_loop.stop()-cdr_loop.start()+1;
	return l;
}

std::string
AntibodyInfo::get_cluster_name(CDRClusterEnum const cluster) {
	return cdr_cluster_manager_->cdr_cluster_enum_to_string(cluster);
}

CDRClusterEnum
AntibodyInfo::get_cluster_enum(const std::string cluster) {
	return cdr_cluster_manager_->cdr_cluster_string_to_enum(cluster);
}

CDRNameEnum
AntibodyInfo::get_cdr_enum_for_cluster(CDRClusterEnum const cluster) {
	std::string cluster_string = cdr_cluster_manager_->cdr_cluster_enum_to_string(cluster);
	utility::vector1< std::string > cluster_vec = utility::string_split(cluster_string, '-');
	return antibody_manager_->cdr_name_string_to_enum(cluster_vec[1]);
}

core::Size
AntibodyInfo::get_cluster_length(CDRClusterEnum const cluster) {
	std::string cluster_string = cdr_cluster_manager_->cdr_cluster_enum_to_string(cluster);
	utility::vector1< std::string > cluster_vec = utility::string_split(cluster_string, '-');
	return utility::string2int(cluster_vec[2]);
}


////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
///				provide fold tree utilities for various purpose                    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

kinematics::FoldTreeCOP AntibodyInfo::setup_simple_fold_tree(
    Size const & jumppoint1,
    Size const & cutpoint,
    Size const & jumppoint2,
    pose::Pose const & pose ) const {
	using namespace kinematics;

	FoldTreeOP f = new FoldTree();
	loops::Loops moveable_region;
	moveable_region.add_loop(loops::Loop(jumppoint1 + 1, jumppoint2 - 1, cutpoint));
	loops::fold_tree_from_loops(pose, moveable_region, *f);
	return f;

}


kinematics::FoldTreeCOP AntibodyInfo::get_FoldTree_AllCDRs (pose::Pose const & pose ) const {
	using namespace kinematics;

	kinematics::FoldTreeOP f = new FoldTree();
	loops::fold_tree_from_loops(pose, *loopsop_having_allcdrs_, *f);
	return f;

} // all_cdr_fold_tree()

///////////////////////////////////////////////////////////////////////////
/// @begin all_cdr_VL_VH_fold_tree
///
/// @brief change to all CDR and VL-VH dock fold tree
///
/// @authors Aroop 07/13/2010
///
/// @last_modified 07/13/2010
///////////////////////////////////////////////////////////////////////////
kinematics::FoldTreeCOP
AntibodyInfo::get_FoldTree_AllCDRs_LHDock( pose::Pose const & pose ) const {

	using namespace kinematics;

	Size nres = pose.total_residue();
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	char second_chain = 'H';
	Size rb_cutpoint(0);

	for ( Size i = 1; i <= nres; ++i ) {
		if( pdb_info->chain( i ) == second_chain) {
			rb_cutpoint = i-1;
			break;
		}
	}

	Size jump_pos1 ( geometry::residue_center_of_mass( pose, 1, rb_cutpoint ) );
	Size jump_pos2 ( geometry::residue_center_of_mass( pose,rb_cutpoint+1, nres ) );
	//TR<<rb_cutpoint<<std::endl;
	//TR<<jump_pos1<<std::endl;
	//TR<<jump_pos2<<std::endl;

	// make sure rb jumps do not reside in the loop region
	for( loops::Loops::const_iterator it= loopsop_having_allcdrs_->begin(), it_end = loopsop_having_allcdrs_->end(); it != it_end; ++it ) {
		if (   jump_pos1 >= ( it->start() - 1 ) && jump_pos1 <= ( it->stop() + 1)   )
			jump_pos1 = it->stop() + 2;
		if (   jump_pos2 >= ( it->start() - 1 ) && jump_pos2 <= ( it->stop() + 1)   )
			jump_pos2 = it->start() - 2;
	}

	// make a simple rigid-body jump first
	FoldTreeOP f = new FoldTree(* setup_simple_fold_tree(jump_pos1,rb_cutpoint,jump_pos2, pose ));

	for( loops::Loops::const_iterator it=loopsop_having_allcdrs_->begin(), it_end=loopsop_having_allcdrs_->end(); it != it_end; ++it ) {
		Size const loop_start ( it->start() );
		Size const loop_stop ( it->stop() );
		Size const loop_cutpoint ( it->cut() );
		Size edge_start(0), edge_stop(0);
		//bool edge_found = false;  // unused ~Labonte
		const FoldTree & f_const = *f;
		Size const num_jump = f_const.num_jump();
		for( FoldTree::const_iterator it2=f_const.begin(), it2_end=f_const.end(); it2 !=it2_end; ++it2 ) {
			//TR<<it2->start()<<std::endl;
			//TR<<it2->stop()<<std::endl;
			edge_start = std::min( it2->start(), it2->stop() );
			edge_stop  = std::max( it2->start(), it2->stop() );
			if ( ! it2->is_jump() && loop_start > edge_start && loop_stop < edge_stop ) {
				//edge_found = true;  // unused ~Labonte
				break;
			}
		}

		f->delete_unordered_edge( edge_start, edge_stop, Edge::PEPTIDE);
		f->add_edge( loop_start-1, loop_stop+1, num_jump+1 );
		f->add_edge( edge_start, loop_start-1, Edge::PEPTIDE );
		f->add_edge( loop_start-1, loop_cutpoint, Edge::PEPTIDE );
		f->add_edge( loop_cutpoint+1, loop_stop+1, Edge::PEPTIDE );
		f->add_edge( loop_stop+1, edge_stop, Edge::PEPTIDE );
	}

	f->reorder(1);
	return f;
}


///////////////////////////////////////////////////////////////////////////
/// @begin LH_A_foldtree
///
/// @brief	Fold tree for snugdock, docks LH with the antigen chains. The function
///			assumes that the coordinates for antigen chains in the input PDB file
///			are right after the antibody heavy chain (which must be named H).The
///			expected order of chains is thus L, H followed by the antigen chains.
///
/// @authors Krishna Praneeth Kilambi 08/14/2012
///
/// @last_modified 08/14/2012
///////////////////////////////////////////////////////////////////////////
kinematics::FoldTree
AntibodyInfo::get_FoldTree_LH_A( pose::Pose const & pose ) const {

	using namespace core;
	using namespace kinematics;

	Size nres = pose.total_residue();
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	char second_chain = 'H';
	Size cutpoint = 0;

	kinematics::FoldTree LH_A_foldtree;

	for ( Size i = 1; i <= nres; ++i ) {
		if(pdb_info->chain(1) != 'L') {
			throw excn::EXCN_Msg_Exception("Chains are not named correctly or are not in the expected order");
			break;
		}
		if( (pdb_info->chain(i) == 'L') && (pdb_info->chain(i) != pdb_info->chain(i+1))) {
			if(pdb_info->chain(i+1) != second_chain) {
				throw excn::EXCN_Msg_Exception("Chains are not named correctly or are not in the expected order");
				break;
			}
		}
		if( (pdb_info->chain(i) == second_chain) && (pdb_info->chain(i) != pdb_info->chain(i+1))) {
			cutpoint = i;
			break;
		}
	}

	Size jump_pos1 ( geometry::residue_center_of_mass( pose, 1, cutpoint ) );
	Size jump_pos2 ( geometry::residue_center_of_mass( pose, cutpoint+1, pose.total_residue() ) );

	//setup fold tree based on cutpoints and jump points
	LH_A_foldtree.clear();
	LH_A_foldtree.simple_tree( pose.total_residue() );
	LH_A_foldtree.new_jump( jump_pos1, jump_pos2, cutpoint);

	Size chain_begin(0), chain_end(0);

	//rebuild jumps between antibody light and heavy chains
	chain_end = cutpoint;
	chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
	while (chain_begin != 1) {
		chain_end = chain_begin-1;
		LH_A_foldtree.new_jump( chain_end, chain_begin, chain_end);
		chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
	}

	//rebuild jumps between all the antigen chains
	chain_begin = cutpoint+1;
	chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	while (chain_end != pose.total_residue()) {
		chain_begin = chain_end+1;
		LH_A_foldtree.new_jump( chain_end, chain_begin, chain_end);
		chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	}

	LH_A_foldtree.reorder( 1 );
	LH_A_foldtree.check_fold_tree();

	return LH_A_foldtree;
}


///////////////////////////////////////////////////////////////////////////
/// @begin L_HA_foldtree
///
/// @brief	Fold tree for LH refinement in snugdock, docks L with H + antigen
///			chains. The function assumes that the coordinates for antigen chains
///			in the input PDB file are right after the antibody heavy chain
///			(which must be named H).The expected order of chains is thus
///			L, H followed by the antigen chains.
///
/// @authors Krishna Praneeth Kilambi 08/14/2012
///
/// @last_modified 08/14/2012
///////////////////////////////////////////////////////////////////////////
kinematics::FoldTree
AntibodyInfo::get_FoldTree_L_HA( pose::Pose const & pose ) const {

	using namespace core;
	using namespace kinematics;

	Size nres = pose.total_residue();
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	char second_chain = 'H';
	Size cutpoint = 0;

	kinematics::FoldTree L_HA_foldtree;

	for ( Size i = 1; i <= nres; ++i ) {
		if(pdb_info->chain(1) != 'L') {
			throw excn::EXCN_Msg_Exception("Chains are not named correctly or are not in the expected order");
			break;
		}
		if( (pdb_info->chain(i) == 'L') && (pdb_info->chain(i) != pdb_info->chain(i+1))) {
			if(pdb_info->chain(i+1) != second_chain) {
				throw excn::EXCN_Msg_Exception("Chains are not named correctly or are not in the expected order");
				break;
			}
		}
		if( (pdb_info->chain(i) == 'L') && (pdb_info->chain(i+1) == second_chain)) {
			cutpoint = i;
			break;
		}
	}

	Size jump_pos1 ( geometry::residue_center_of_mass( pose, 1, cutpoint ) );
	Size jump_pos2 ( geometry::residue_center_of_mass( pose, cutpoint+1, pose.conformation().chain_end( pose.chain(cutpoint+1) ) ) );

	//setup fold tree based on cutpoints and jump points
	L_HA_foldtree.clear();
	L_HA_foldtree.simple_tree( pose.total_residue() );
	L_HA_foldtree.new_jump( jump_pos1, jump_pos2, cutpoint);

	Size chain_begin(0), chain_end(0);

	//rebuild jumps between heavy chain and antigen chains
	chain_begin = cutpoint+1;
	chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	while (chain_end != pose.total_residue()) {
		chain_begin = chain_end+1;
		L_HA_foldtree.new_jump( chain_end, chain_begin, chain_end);
		chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	}

	L_HA_foldtree.reorder( 1 );
	L_HA_foldtree.check_fold_tree();

	return L_HA_foldtree;

}

///////////////////////////////////////////////////////////////////////////
/// @begin LA_H_foldtree
///
/// @brief	Fold tree for LH refinement in snugdock, docks L + antigen chains
///			with H. The function assumes that the coordinates for antigen chains
///			in the input PDB file are right after the antibody heavy chain
///			(which must be named H).The expected order of chains is thus
///			L, H followed by the antigen chains.
///
/// @authors Krishna Praneeth Kilambi 08/14/2012
///
/// @last_modified 08/14/2012
///////////////////////////////////////////////////////////////////////////
kinematics::FoldTree
AntibodyInfo::get_FoldTree_LA_H( pose::Pose const & pose ) const {

	using namespace core;
	using namespace kinematics;

	Size nres = pose.total_residue();
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	char second_chain = 'H';
	Size cutpoint = 0;
	bool lchain_jump = false;

	kinematics::FoldTree LA_H_foldtree ;

	for ( Size i = 1; i <= nres; ++i ) {
		if(pdb_info->chain(1) != 'L') {
			throw excn::EXCN_Msg_Exception("Chains are not named correctly or are not in the expected order");
			break;
		}
		if( (pdb_info->chain(i) == 'L') && (pdb_info->chain(i) != pdb_info->chain(i+1))) {
			if(pdb_info->chain(i+1) != second_chain) {
				throw excn::EXCN_Msg_Exception("Chains are not named correctly or are not in the expected order");
				break;
			}
		}
		if( (pdb_info->chain(i) == 'L') && (pdb_info->chain(i+1) == second_chain)) {
			cutpoint = i;
			break;
		}
	}

	Size jump_pos1 ( geometry::residue_center_of_mass( pose, 1, cutpoint ) );
	Size jump_pos2 ( geometry::residue_center_of_mass( pose, cutpoint+1, pose.conformation().chain_end( pose.chain(cutpoint+1) ) ) );

	//setup fold tree based on cutpoints and jump points
	LA_H_foldtree.clear();
	LA_H_foldtree.simple_tree( pose.total_residue() );
	LA_H_foldtree.new_jump( jump_pos1, jump_pos2, cutpoint);

	Size chain_begin(0), chain_end(0);

	//rebuild jumps between the light chain and antigen chains
	chain_begin = cutpoint+1;
	chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	while (chain_end != pose.total_residue()) {
		chain_begin = chain_end+1;
		if (!lchain_jump) {
			LA_H_foldtree.new_jump( pose.conformation().chain_end( pose.chain(1) ), chain_begin, chain_end);
			lchain_jump = true;
		} else {
			LA_H_foldtree.new_jump( chain_end, chain_begin, chain_end);
		}
		chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
	}

	LA_H_foldtree.reorder( 1 );
	LA_H_foldtree.check_fold_tree();

	return LA_H_foldtree;


}

kinematics::MoveMap
AntibodyInfo::get_MoveMap_for_Loops(pose::Pose const & pose,
                                    loops::Loops const & the_loops,
                                    bool const & bb_only,
                                    bool const & include_nb_sc,
                                    Real const & nb_dist) const {
	kinematics::MoveMap move_map ;

	move_map.clear();
	move_map.set_chi( false );
	move_map.set_bb( false );
	utility::vector1< bool> bb_is_flexible( pose.total_residue(), false );
	utility::vector1< bool> sc_is_flexible( pose.total_residue(), false );

	select_loop_residues( pose, the_loops, false/*include_neighbors*/, bb_is_flexible, nb_dist);
	move_map.set_bb( bb_is_flexible );
	if (bb_only==false) {
		select_loop_residues( pose, the_loops, include_nb_sc/*include_neighbors*/, sc_is_flexible, nb_dist);
		move_map.set_chi( sc_is_flexible );
	}
	for( Size ii = 1; ii <= the_loops.num_loop(); ++ii ) {
		move_map.set_jump( ii, false );
	}

	return move_map;
}

kinematics::MoveMap
AntibodyInfo::get_MoveMap_for_AllCDRsSideChains_and_H3backbone(pose::Pose const & pose,
        bool const & include_nb_sc ,
        Real const & nb_dist ) const {

	kinematics::MoveMap move_map ;

	move_map.clear();
	move_map.set_chi( false );
	move_map.set_bb( false );
	utility::vector1< bool> bb_is_flexible( pose.total_residue(), false );
	utility::vector1< bool> sc_is_flexible( pose.total_residue(), false );

	TR<<"start: "<<get_CDR_start(h3, pose)<<std::endl;
	TR<<"end  : "<<get_CDR_end(h3, pose)<<std::endl;
	for(Size ii=get_CDR_start(h3, pose); ii<=get_CDR_end(h3, pose); ii++) {
		bb_is_flexible[ii] = true;
		TR<<"Setting Residue "<< ii<<" to be true"<<std::endl;
	}
	move_map.set_bb( bb_is_flexible );

	select_loop_residues( pose, *(get_AllCDRs_in_loopsop()), include_nb_sc, sc_is_flexible, nb_dist);
	move_map.set_chi( sc_is_flexible );

	for( Size ii = 1; ii <= (get_AllCDRs_in_loopsop())->num_loop(); ++ii ) {
		move_map.set_jump( ii, false );
	}


	return move_map;
}

void
AntibodyInfo::add_CDR_to_MoveMap(pose::Pose const & pose,
                                 kinematics::MoveMapOP movemap,
                                 CDRNameEnum const & cdr_name,
                                 bool const & bb_only,
                                 bool const & include_nb_sc,
                                 Real const & nb_dist) const {

	core::Size start = get_CDR_start(cdr_name, pose);
	core::Size stop =  get_CDR_end(cdr_name, pose);
	core::Size cut = (stop-start+1)/2+start;

	loops::Loop cdr_loop = loops::Loop(start, stop, cut);
	utility::vector1< bool> bb_is_flexible( pose.total_residue(), false );
	select_loop_residues( pose, cdr_loop, include_nb_sc/*include_neighbors*/, bb_is_flexible, nb_dist);
	for (core::Size i=1; i<=pose.total_residue(); ++i) {
		if (bb_is_flexible[i]) {
			movemap->set_bb(i, true);
			if (! bb_only) {
				movemap->set_chi(i, true);
			}
		}
	}
}

kinematics::MoveMap
AntibodyInfo::get_MoveMap_for_LoopsandDock(pose::Pose const & pose,
        loops::Loops const & the_loops,
        bool const & bb_only,
        bool const & include_nb_sc,
        Real const & nb_dist) const {

	kinematics::MoveMap move_map = get_MoveMap_for_Loops(pose, the_loops, bb_only, include_nb_sc, nb_dist);

	move_map.set_jump( 1, true );
	for( Size ii = 2; ii <= the_loops.num_loop() +1 ; ++ii ) {
		move_map.set_jump( ii, false );
	}

	return move_map;
}




//JQX: doesn't matter only antibody or antibody-antigen complex, just include CDRs and their neighbors
pack::task::TaskFactoryOP
AntibodyInfo::get_TaskFactory_AllCDRs(pose::Pose & pose)  const {

	vector1< bool> sc_is_packable( pose.total_residue(), false );
	select_loop_residues( pose, *loopsop_having_allcdrs_, true/*include_neighbors*/, sc_is_packable);

	using namespace pack::task;
	using namespace pack::task::operation;
	// selecting movable c-terminal residues
	ObjexxFCL::FArray1D_bool loop_residues( pose.total_residue(), false );
	for( Size i = 1; i <= pose.total_residue(); ++i ) {
		loop_residues(i) = sc_is_packable[i];
	} // check mapping

	using namespace protocols::toolbox::task_operations;
	pack::task::TaskFactoryOP tf ;
	tf= setup_packer_task(pose);
	//	tf->push_back( new RestrictToInterface(loop_residues) ); //JQX: not sure why we use loop_residues, in stead of sc_is_packable
	tf->push_back( new RestrictToInterface(sc_is_packable) );



	//pack::task::PackerTaskOP my_task2(tf->create_task_and_apply_taskoperations(pose));
	//TR<<*my_task2<<std::endl; //exit(-1);

	return tf;
}

pack::task::TaskFactoryOP
AntibodyInfo::get_TaskFactory_OneCDR(pose::Pose & pose, CDRNameEnum const & cdr_name)  const {
	vector1< bool> sc_is_packable( pose.total_residue(), false );

	select_loop_residues( pose, *get_CDR_in_loopsop(cdr_name), true/*include_neighbors*/, sc_is_packable);
	using namespace protocols::toolbox::task_operations;

	pack::task::TaskFactoryOP tf ;
	tf= setup_packer_task(pose);
	tf->push_back( new RestrictToInterface(sc_is_packable) );

	return tf;
}

/// TODO:
//JQX: make Daisuke's code compatible with my code
//JAB: Expand using HMMER antibody/cdr identification

///////////////////////////////////////////////////////////////////////////
/// @author: Daisuke Kuroda (dkuroda1981@gmail.com) 06/18/2012
///
/// @brief: Identify 3 CDRs from a sequence. Automatically judge heavy or light chains (I hope!).
///         The input can be either a light chain, a heavy chain or another sequence.
///
/// @last_modified 08/28/2012 by DK
///////////////////////////////////////////////////////////////////////////
void
AntibodyInfo::identify_CDR_from_a_sequence(std::string const & querychain) {

	int l1found = 0, l2found = 0, l3found = 1, h1found = 1, h2found = 0, h3found = 1; // 0 if exist; otherwise 1.
	int lenl1 = 0, lenl2 = 0, lenl3 = 0, lenh1 = 0, lenh2 = 0, lenh3 = 0;
	int posl1_s = 0, posl1_e = 0, posl2_s = 0, posl2_e = 0, posl3_s = 0, posl3_e = 0;
	int posh1_s = 0, posh1_e = 0, posh2_s = 0, posh2_e = 0, posh3_s = 0, posh3_e = 0;
	// int i = 0, // Unused variable causes warning.
	int k = 0, l = 0, m = 0, n = 0;

	int pos_fr1_s = 0, pos_fr1_e = 0, pos_fr2_s = 0, pos_fr2_e = 0;
	int pos_fr3_s = 0, pos_fr3_e = 0, pos_fr4_s = 0, pos_fr4_e = 0;
	int len_fr1 = 0, len_fr2 = 0, len_fr3 = 0, len_fr4 = 0;

	int len;
	std::string check;

	std::string seql1, seql2, seql3, seqh1, seqh2, seqh3;
	std::string frl3, frh1, frh3;

	std::string seq_fr1, seq_fr2, seq_fr3, seq_fr4;

	// For L3:	[FVI]-[GAEV]-X-[GY]
	std::string p1_l3[] = {"F","V","I"};
	std::string p2_l3[] = {"G","A","E","V"};
	std::string p3_l3[] = {"G","A","P","C","D","E","Q","N","R","K","H","W","Y","F","M","T","V","I","S","L"};
	std::string p4_l3[] = {"G","Y"};

	// For H1
	std::string p1_h1[] = {"W","L"};
	std::string p2_h1[] = {"I","V","F","Y","A","M","L","N","G","E","W"};
	std::string p3_h1[] = {"R","K","Q","V","N","C"};
	std::string p4_h1[] = {"Q","K","H","E","L","R"};

	// For H3	W-[GA]-X-[DRG]
	// std::string p1_h3[] = {"W","V"}; // Unused variable causes warning.
	std::string p2_h3[] = {"G","A","C"};
	std::string p3_h3[] = {"A","P","C","D","E","Q","N","R","K","H","W","Y","F","M","T","V","I","S","L","G"};
	std::string p4_h3[] = {"G","R","D","Q","V"};

	// Input sequence is here
	std::string querychain2(querychain, 0,140);
	std::string querychain3(querychain, 0,110);

	std::string querychain_first(querychain, 0, 50);
	std::string querychain_last(querychain, 70);

	//TR << querychain2 << endl;
	//TR << querychain_last << endl;

	len = querychain.length();
	if(len < 130) {
		check = "Fv";
	} else if(len < 250) {
		check = "Fab";
	} else {
		check = "Weird";
	}

	//TR << "*** Query sequence ***" << endl;
	//TR << querychain << endl;
	//TR << endl;

	/*****************************************************/
	/***************** Is it light chain? ****************/
	/*****************************************************/
	/* L1 search Start */
	if(querychain_first.find("WYL") != std::string::npos) {
		posl1_e = querychain_first.find("WYL") - 1;
	} else if(querychain_first.find("WLQ") != std::string::npos) {
		posl1_e = querychain_first.find("WLQ") - 1;
	} else if(querychain_first.find("WFQ") != std::string::npos) {
		posl1_e = querychain_first.find("WFQ") - 1;
	} else if(querychain_first.find("WYQ") != std::string::npos) {
		posl1_e = querychain_first.find("WYQ") - 1;
	} else if(querychain_first.find("WYH") != std::string::npos) {
		posl1_e = querychain_first.find("WYH") - 1;
	} else if(querychain_first.find("WVQ") != std::string::npos) {
		posl1_e = querychain_first.find("WVQ") - 1;
	} else if(querychain_first.find("WVR") != std::string::npos) {
		posl1_e = querychain_first.find("WVR") - 1;
	} else if(querychain_first.find("WWQ") != std::string::npos) {
		posl1_e = querychain_first.find("WWQ") - 1;
	} else if(querychain_first.find("WVK") != std::string::npos) {
		posl1_e = querychain_first.find("WVK") - 1;
	} else if(querychain_first.find("WLL") != std::string::npos) {
		posl1_e = querychain_first.find("WLL") - 1;
	} else if(querychain_first.find("WFL") != std::string::npos) {
		posl1_e = querychain_first.find("WFL") - 1;
	} else if(querychain_first.find("WVF") != std::string::npos) {
		posl1_e = querychain_first.find("WVF") - 1;
	} else if(querychain_first.find("WIQ") != std::string::npos) {
		posl1_e = querychain_first.find("WIQ") - 1;
	} else if(querychain_first.find("WYR") != std::string::npos) {
		posl1_e = querychain_first.find("WYR") - 1;
	} else if(querychain_first.find("WNQ") != std::string::npos) {
		posl1_e = querychain_first.find("WNQ") - 1;
	} else if(querychain_first.find("WHL") != std::string::npos) {
		posl1_e = querychain_first.find("WHL") - 1;
	} else if(querychain_first.find("WYM") != std::string::npos) {	// Add 06/10/2012
		posl1_e = querychain_first.find("WYM") - 1;
	} else {
		l1found = 1;
	}

	if(l1found != 1) {
		posl1_s   = querychain_first.find("C") + 1;
		lenl1     = posl1_e - posl1_s + 1;
		seql1     = querychain_first.substr(posl1_s,lenl1);

		pos_fr1_s = 0;
		pos_fr1_e = posl1_s - 1;
		len_fr1   = pos_fr1_e - pos_fr1_s + 1;
		seq_fr1   = querychain_first.substr(pos_fr1_s,len_fr1);
	}
	/* L1 search Finish */

	/* L2 search start */
	if(l1found != 1) {
		posl2_s   = posl1_e + 16;
		posl2_e   = posl2_s + 6;
		lenl2     = posl2_e - posl2_s + 1;
		seql2     = querychain.substr(posl2_s,lenl2);

		pos_fr2_s = posl1_e + 1;
		pos_fr2_e = posl2_s - 1;
		len_fr2   = pos_fr2_e - pos_fr2_s + 1;
		seq_fr2   = querychain.substr(pos_fr2_s,len_fr2);
	} else {
		l2found = 1;
	}
	/* L2 search end */

	/* L3 search Start */
	//string p1_l3[] = {"F","V","I"};
	//string p2_l3[] = {"G","A","E","V"};
	//string p3_l3[] = {"G","A","P","C","D","E","Q","N","R","K","H","W","Y","F","M","T","V","I","S","L"};
	//string p4_l3[] = {"G","Y"};
	for(l = 0; l < 3; ++l) {
		for(m = 0; m < 4; ++m) {
			for(n = 0; n < 2; ++n) {
				for(k = 0; k < 20; ++k) {
					//frl3 = "FG" + p3_l3[k] + "G";
					//frl3 = p1_l3[l] + "G" + p3_l3[k] + "G";			//[VF]GXG
					frl3 = p1_l3[l] + p2_l3[m] + p3_l3[k] + p4_l3[n];	//[VF][AG]XG

					if(querychain3.find(frl3, 80) != std::string::npos) {
						posl3_e = querychain3.find(frl3,80) - 1;
						posl3_s = querychain3.find("C",80) + 1;
						lenl3   = posl3_e - posl3_s + 1;

						//TR << frl3 << "\t" << posl3_s << "\t" << lenl3 << endl;

						seql3   = querychain3.substr(posl3_s,lenl3);

						if(seql3.length() > 4) {
							l3found = 0;
							break;
						} else {
							l3found = 1;
						}

						pos_fr3_s = posl2_e + 1;
						pos_fr3_e = posl3_s - 1;
						pos_fr4_s = posl3_e + 1;
						pos_fr4_e = pos_fr4_s + 5;
						len_fr3   = pos_fr3_e - pos_fr3_s + 1;
						len_fr4   = pos_fr4_e - pos_fr4_s + 1;
						seq_fr3   = querychain3.substr(pos_fr3_s,len_fr3);
						seq_fr4   = querychain3.substr(pos_fr4_s,len_fr4);
					}
				}

				l = 3;
				m = 4;
				n = 2;
				k = 20;
			}
		}
	}
	/* L3 search Finish */

	/*****************************************************/
	/***************** Is it heavy chain? ****************/
	/*****************************************************/
	if(l1found == 1 || l2found == 1 || l3found == 1) {
		/* H1 search Start */
		//string p1_h1[] = {"W","L"};
		//string p2_h1[] = {"I","V","F","Y","A","M","L","N","G","E"};
		//string p3_h1[] = {"R","K","Q","V","N","C"};
		//string p4_h1[] = {"Q","K","H","E","L","R"};
		for(n = 0; n < 2; ++n) {
			for(l = 0; l < 6; ++l) {
				for(m = 0; m < 6; ++m) {
					for(k = 0; k < 10; ++k) {
						//frh1 = "W" + p2_h1[k] + p3_h1[l] + p4_h1[m];
						frh1 = p1_h1[n] + p2_h1[k] + p3_h1[l] + p4_h1[m];

						if(querychain_first.find(frh1, 0) != std::string::npos) {
							posh1_e = querychain_first.find(frh1, 0) - 1;
							h1found = 0;
							n = 2;
							l = 6;
							m = 6;
							k = 10;
						}
					}
				}
			}
		}

		if(h1found != 1) {
			posh1_s  = querychain_first.find("C") + 4;
			lenh1    = posh1_e - posh1_s + 1;
			seqh1    = querychain_first.substr(posh1_s, lenh1);

			pos_fr1_s = 0;
			pos_fr1_e = posh1_s - 1;
			len_fr1   = pos_fr1_e - pos_fr1_s + 1;
			seq_fr1   = querychain_first.substr(pos_fr1_s,len_fr1);
		}
		/* H1 search Finish */

		/* H3 search Start */
		//string p1_h3[] = {"W","V"};
		//string p2_h3[] = {"G","A","C"};
		//string p3_h3[] = {"A","P","C","D","E","Q","N","R","K","H","W","Y","F","M","T","V","I","S","L","G"};
		//string p4_h3[] = {"G","R","D","Q","V"};
		for(m = 0; m < 3; ++m) {
			for(l = 0; l < 5; ++l) {
				for(k = 0; k < 20; ++k) {
					//frh3 = "WG" + p3_h3[k] + p4_h3[l];
					frh3 = "W" + p2_h3[m] + p3_h3[k] + p4_h3[l];

					if(querychain2.find(frh3, 80) != std::string::npos) {
						posh3_e = querychain2.find(frh3,80) - 1;
						h3found = 0;
						m = 3;
						l = 5;
						k = 20;
					}
				}
			}
		}

		if(querychain2.find("C", 80) != std::string::npos) {
			posh3_s = querychain2.find("C", 80) + 3;
		} else {
			h3found = 1;
		}

		if(h3found != 1) {
			lenh3 = posh3_e - posh3_s + 1;
			seqh3 = querychain2.substr(posh3_s,lenh3);
		}
		/* H3 search Finish */

		/* H2 search start */
		if(h1found != 1 && h3found != 1) {
			posh2_s = posh1_e + 15;
			posh2_e = posh3_s - 33;
			lenh2   = posh2_e - posh2_s + 1;
			seqh2   = querychain.substr(posh2_s,lenh2);

			pos_fr2_s = posh1_e + 1;
			pos_fr2_e = posh2_s - 1;
			pos_fr3_s = posh2_e + 1;
			pos_fr3_e = posh3_s - 1;
			pos_fr4_s = posh3_e + 1;
			pos_fr4_e = pos_fr4_s + 5;
			len_fr2   = pos_fr2_e - pos_fr2_s + 1;
			len_fr3   = pos_fr3_e - pos_fr3_s + 1;
			len_fr4   = pos_fr4_e - pos_fr4_s + 1;
			seq_fr2   = querychain.substr(pos_fr2_s,len_fr2);
			seq_fr3   = querychain.substr(pos_fr3_s,len_fr3);
			seq_fr4   = querychain.substr(pos_fr4_s,len_fr4);
		}
		/* H2 search end */
	}

	if(l1found == 0 && l2found == 0 && l3found == 0) {
		TR << lenl1 << "\t" << lenl2 << "\t" << lenl3 << "\t";
		TR << seql1 << "\t" << seql2 << "\t" << seql3 << "\t" << check << "\tLIGHT" << std::endl;
	} else if(h1found == 0 && h2found == 0 && h3found == 0) {
		TR << lenh1 << "\t" << lenh2 << "\t" << lenh3 << "\t";
		TR << seqh1 << "\t" << seqh2 << "\t" << seqh3 << "\t" << check << "\tHEAVY" << std::endl;
	} else if(l1found == 0 && l2found == 0 && l3found == 0 && h1found == 0 && h2found == 0 && h3found == 0) {
		TR << lenl1 << "\t" << lenl2 << "\t" << lenl3 << "\t";
		TR << lenh1 << "\t" << lenh2 << "\t" << lenh3 << "\t";
		TR << seql1 << "\t" << seql2 << "\t" << seql3 << "\t" << check << "\tLIGHT" << "\t";
		TR << seqh1 << "\t" << seqh2 << "\t" << seqh3 << "\t" << check << "\tHEAVY" << std::endl;
	} else {
		TR << "Some CDRs seem to be missing!\t" << querychain << "\t";
		TR << lenh1 << "\t" << lenh2 << "\t" << lenh3 << "\t";
		TR << "H1:" << seqh1 << "\tH2: " << seqh2 << "\tH3 " << seqh3 << "\t" << posh3_s << std::endl;
	}

	// return 0;
}




vector1<char> AntibodyInfo::get_CDR_sequence_with_stem(CDRNameEnum const & cdr_name,
        Size left_stem ,
        Size right_stem ) const {
	vector1<char> sequence;
	loops::Loop the_loop = get_CDR_loop(cdr_name);

	for (Size i=the_loop.start()-left_stem; i<=the_loop.stop()+right_stem; ++i) {
		std::map< Size, char >::const_iterator iter(sequence_map_.find(i)); //To keep const
		sequence.push_back(iter->second);
	}
	return sequence;
}

scoring::ScoreFunctionCOP get_Pack_ScoreFxn(void) {
	static scoring::ScoreFunctionOP pack_scorefxn = 0;
	if(pack_scorefxn == 0) {
		pack_scorefxn = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}
	return pack_scorefxn;
}
scoring::ScoreFunctionCOP get_Dock_ScoreFxn(void) {
	static scoring::ScoreFunctionOP dock_scorefxn = 0;
	if(dock_scorefxn == 0) {
		dock_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
		dock_scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
		dock_scorefxn->set_weight( core::scoring::overlap_chainbreak, 10./3. );
	}
	return dock_scorefxn;
}
scoring::ScoreFunctionCOP get_LoopCentral_ScoreFxn(void) {
	static scoring::ScoreFunctionOP loopcentral_scorefxn = 0;
	if(loopcentral_scorefxn == 0) {
		loopcentral_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
		loopcentral_scorefxn->set_weight( scoring::chainbreak, 10./3. );
	}
	return loopcentral_scorefxn;
}
scoring::ScoreFunctionCOP get_LoopHighRes_ScoreFxn(void) {
	static scoring::ScoreFunctionOP loophighres_scorefxn = 0;
	if(loophighres_scorefxn == 0) {
		loophighres_scorefxn = scoring::getScoreFunction();
		loophighres_scorefxn->set_weight( scoring::chainbreak, 1.0 );
		loophighres_scorefxn->set_weight( scoring::overlap_chainbreak, 10./3. );
		loophighres_scorefxn->set_weight(scoring::dihedral_constraint, 1.0);
	}
	return loophighres_scorefxn;
}

AntibodyEnumManagerOP &
AntibodyInfo::get_antibody_enum_manager() {
	return antibody_manager_;
}

CDRClusterEnumManagerOP &
AntibodyInfo::get_cdr_cluster_enum_manager() {
	return cdr_cluster_manager_;
}


/// @details  Show the complete setup of the docking protocol
void AntibodyInfo::show( std::ostream & out ) {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const AntibodyInfo & ab_info )  {

	using namespace ObjexxFCL::fmt;
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << A( 47, "Rosetta Antibody Info" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;

	out << line_marker << "             Antibody Type:";
	if(ab_info.is_camelid()) {
		out << "  Camelid Antibody"<< std::endl;
	} else                    {
		out << "  Regular Antibody"<< std::endl;
	}

	out << line_marker << " Predict H3 Cterminus Base:";
	out <<"  "<<ab_info.antibody_manager_->h3_base_type_enum_to_string(ab_info.predicted_H3_base_type_)<<std::endl;

	out << line_marker << space( 74 ) << std::endl;
	for (CDRNameEnum i=start_cdr_loop; i<=ab_info.total_cdr_loops_; i=CDRNameEnum(i+1) ) {
		out << line_marker << " "+ab_info.get_CDR_Name(i)+" info: "<<std::endl;
		out << line_marker << "           length:  "<< ab_info.get_CDR_loop(i).length() <<std::endl;
		out << line_marker << "         sequence:  ";
		for (Size j=1; j<=ab_info.get_CDR_sequence_with_stem(i,0,0).size(); ++j) {
			out << ab_info.get_CDR_sequence_with_stem(i,0,0)[j] ;
		}
		out <<std::endl;
		out << line_marker << "        loop_info:  "<< ab_info.get_CDR_loop(i)<<std::endl;
	}
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}



} // namespace antibody
} // namespace protocols
