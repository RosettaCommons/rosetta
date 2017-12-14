// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/evolution/AlignmentCleanerTools.hh
/// @author Christoffer Norn (ch.norn@gmail.com)

#include <string>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>


namespace protocols {
namespace evolution {
namespace AlignmentCleanerTools {

void thread_sequence_on_pose( core::pose::Pose & pose, std::string const & thread_seq, core::scoring::ScoreFunctionOP scorefxn );

std::tuple< std::string, std::string > indel_motif(std::string const & aln_seq, core::Size const motif_radius, core::Size const aln_resi, std::string const & pose_ss_aln);

core::Real indel_motif_seq_id( std::string const & motif1, std::string const & motif2 );

std::string short_ss_loop_filter( std::string ss, core::Size min_loop_length);

} // AlignmentCleanerTools
} // evolution
} // protocols
