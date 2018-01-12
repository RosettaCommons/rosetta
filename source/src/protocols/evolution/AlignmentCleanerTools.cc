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

#include <protocols/evolution/AlignmentCleanerTools.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>



namespace protocols {
namespace evolution {
namespace AlignmentCleanerTools {

static basic::Tracer TR( "protocols.evolution.AlignmentCleanerTools" );


void
thread_sequence_on_pose(core::pose::Pose & pose, std::string const & thread_seq, core::scoring::ScoreFunctionOP scorefxn)
{
	// Now we thread on the sequence
	using namespace protocols::toolbox::task_operations;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	ThreadSequenceOperationOP tso( new ThreadSequenceOperation );
	tso->target_sequence( thread_seq );
	tso->allow_design_around(false);

	TaskFactoryOP tf;
	tf = TaskFactoryOP( new TaskFactory );
	tf->push_back(tso);
	PackerTaskOP ptask = tf->create_task_and_apply_taskoperations(pose);

	protocols::minimization_packing::PackRotamersMoverOP pack;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		utility_exit_with_message("Not developed for symmetry!");
	} else {
		pack = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover( scorefxn, ptask ) );
	}
	pack->apply( pose );
	(*scorefxn)(pose);
}


std::tuple< std::string, std::string >
indel_motif(std::string const & aln_seq, core::Size const motif_radius, core::Size const aln_resi, std::string const & pose_ss_aln)
{
	// Find the motif sequence downstream
	core::Size first_pos = aln_resi;
	core::Size non_gap_loop_counter = 0;
	for ( core::Size i=aln_resi; i + 1 > 0; --i, --first_pos ) {
		if ( first_pos == 0 ) break;
		if ( aln_seq[i] != '-' && pose_ss_aln[i] != 'L' ) ++non_gap_loop_counter;
		if ( non_gap_loop_counter > motif_radius ) break; // we want it to be > as we start at position
	}

	// Find the motif sequence upstream
	core::Size last_pos = aln_resi;
	non_gap_loop_counter = 0;
	for ( core::Size i = aln_resi; i < aln_seq.length(); ++i, ++last_pos ) {
		if ( aln_seq[i] != '-' && pose_ss_aln[i] != 'L' ) ++non_gap_loop_counter;
		if ( non_gap_loop_counter > motif_radius ) break; // we want it to be > as we start at position
	}

	// Convert to binary motif
	core::Size const len = last_pos - first_pos + 1;
	std::string const aa_motif = aln_seq.substr(first_pos, len);
	std::string binary_motif;
	for ( const char& c: aa_motif ) {
		if ( c == '-' ) binary_motif += '-';
		else binary_motif += 'X';
	}
	return std::make_tuple(aa_motif, binary_motif);
}

core::Real
indel_motif_seq_id( std::string const & motif1, std::string const & motif2)
{
	core::Size len = motif1.length();
	core::Size identities = 0;
	for ( core::Size i=0; i < len; ++i ) {
		if ( motif1[i] == motif2[i] ) ++identities;
	}
	core::Real seq_id = identities /(core::Real) len;
	return seq_id;
}

std::string
short_ss_loop_filter( std::string ss, core::Size min_loop_length)
{
	// Input
	//          LLHHLLHLLLHHH, min_loop_length = 2
	// output
	//          RRHHRRHLLLHHH

	// Check input
	core::Size const window_size = min_loop_length + 2;
	if ( window_size > ss.size() ) {
		utility_exit_with_message("For filtering the secondary structure string, the sliding window (min_loop_length+2) must be smaller than the pose length");
	}

	// Mutate ss where loop is short (<=min_loop_dist) to R
	std::string const min_loop(min_loop_length, min_loop_length);
	for ( core::Size i=0; i < ss.size() - window_size + 1; ++i ) {
		std::string window_ss = ss.substr(i, window_size);

		if ( i == 0 && window_ss[0] == 'L' && window_ss[window_size-2] != 'L' ) { // if Nterm-LLHH...?
			for ( core::Size j = 0; j < window_size - 2; ++j ) {
				if ( ss[j] == 'L' ) ss[j] = 'R';
			}
		} else if ( i == ss.size() - window_size /*=last window*/ && window_ss[window_size-1] == 'L' && window_ss[1] != 'L' ) { // if ...LHLL-cterm?
			for ( core::Size j = ss.size() - window_size + 2; j < ss.size(); ++j ) {
				if ( ss[j] == 'L' ) ss[j] = 'R';
			}
		} else if ( window_ss[0] != 'L' && window_ss[window_size-1] != 'L' ) { // if mid sequence: ...HLLH...
			for ( core::Size j=i+1; j < i + window_size - 1; ++j ) {
				if ( ss[j] == 'L' ) ss[j] = 'R';
			}
		}
	}

	return ss;
}

} // AlignmentCleanerTools
} // evolution
} // protocols
