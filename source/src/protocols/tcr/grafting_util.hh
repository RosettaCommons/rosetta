// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/grafting_util.hh
/// @brief grafting utility functions for protocols/tcr
/// @author Ragul (ragul@umd.edu)

#ifndef INCLUDED_protocols_tcr_grafting_util_hh
#define INCLUDED_protocols_tcr_grafting_util_hh

#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <protocols/tcr/util.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/Pose.hh>


namespace protocols {
namespace tcr {

/// @brief grafting of cdr segments of template structure
/// @details for both alpha and beta chain cdr segment
void graft_cdr(TCRmodel::tmpltinfo &currtmplt, std::string const &template_pdb_path, core::Size const &cdr_start, core::Size const &cdr_end);

/// @brief grafting of Framework/germline template structure
/// @details for both alpha and beta chain framework and germline segment
void graft_framework(TCRmodel::tmpltinfo &currtmplt, std::string const &template_pdb_path, TCRseqInfo::tcrposi const &aho_pos);

/// @brief combine graft cdr structures and framework structure to build crude model
/// @details uses CCDEndsGraftMover with default overhang residue size of 3
void graft_cdr_to_fw(core::pose::Pose &scafold, core::pose::Pose const &cdr_piece, core::Size const &start, core::Size const &end, core::Size &nter_overhang, core::Size &cter_overhang );

} //protocols
} //tcr


#endif //protocols/tcr_grafting_util_hh

