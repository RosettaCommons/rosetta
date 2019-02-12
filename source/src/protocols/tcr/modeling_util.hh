// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/grafting_util.hh
/// @brief modeling utility functions for protocols/tcr
/// @author Ragul (ragul@umd.edu)

#ifndef INCLUDED_protocols_tcr_modeling_util_hh
#define INCLUDED_protocols_tcr_modeling_util_hh

#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRloopRefine.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <protocols/tcr/util.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/loops_main.hh>

namespace protocols {
namespace tcr {

/// @brief orient TCR alpha/beta chains based on the orientation templates
void orient_tcr_chain( core::pose::Pose &fr_piece, core::pose::Pose const &ori_tmplt_chain, TCRseqInfo::tcrposi const &selepos);

/// @brief match template seq with target seq during modeling
/// @details uses protocols::simple_moves::MutateResidue to mutate residues
void match_template_and_target_sequence( core::pose::Pose const &scafold, std::string const &targetseq, core::pose::Pose &posecopy );

/// @brief loop rebuild function
/// @details rebuild loop segment of input pose using given loop start and end position
/// @details takes options loop remodel mover type
void ind_remodel_tcr_cdr_loops( core::pose::Pose &tcr_pose, core::Size const loopstart, core::Size const loopend );

/// @brief refine all cdr loops
/// @details refinement all cdr loops CDR1, CDR2, CDR3 and HV4 loops
void refine_all_cdr_loops ( core::pose::Pose &tcr_pose, TCRseqInfo::tcrposi const &a_posi, TCRseqInfo::tcrposi const &b_posi);

/// @brief loop remodel function for CDR loop
void remodel_tcr_cdr_loops ( core::pose::Pose &tcr_pose, protocols::loops::LoopsOP const &cdr3loops, bool refine );

/// @brief loop refine function for CDR loop
void refine_tcr_cdr_loops ( core::pose::Pose &tcr_pose, protocols::loops::LoopsOP const &cdrloops );

/// @brief loop refine function for CDR3 alpha and CDR3 beta loops
void refine_cdr3_loops ( core::pose::Pose &tcr_pose, antibody::grafting::CDR_Bounds const &acdr3, antibody::grafting::CDR_Bounds const &bcdr3);

/// @brief loop refine function for CDR alpha/beta loop
void refine_cdr_loop ( core::pose::Pose &tcr_pose, antibody::grafting::CDR_Bounds const &cdr);

} //protocols
} //tcr


#endif //protocols/tcr_grafting_util_hh

