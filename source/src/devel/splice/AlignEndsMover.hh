// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AlignEndsMover.hh

#ifndef INCLUDED_devel_splice_AlignEndsMover_hh
#define INCLUDED_devel_splice_AlignEndsMover_hh

#include <devel/splice/AlignEndsMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <algorithm>
#include <core/kinematics/Jump.hh>

// C++ Headers
namespace devel {
namespace splice {

class AlignEndsMover : public protocols::moves::Mover {
public:
	AlignEndsMover();
	~AlignEndsMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	core::Real distance_threshold() const{ return distance_threshold_; }
	void distance_threshold( core::Real const r ){ distance_threshold_ = r;}

	core::Size neighbors() const{ return neighbors_; }
	void neighbors( core::Size n ){ neighbors_ = n; }

	core::Size N_terminal_count() const{ return N_terminal_count_; }
	void N_terminal_count( core::Size n ){ N_terminal_count_ = n; }

	bool odd() const{ return odd_; }
	void odd( bool const b ) { odd_ = b; }

	bool even() const{ return even_; }
	void even( bool const b ) { even_ = b; }

	core::pose::PoseOP template_pose() const;
	void template_pose( core::pose::PoseOP p );

	core::Size stagger() const{ return stagger_; }
	void stagger( core::Size const s ){ stagger_ = s; }

	core::Size strand_length() const{ return strand_length_; }
	void strand_length( core::Size const s ){ strand_length_ = s; }

	core::Size max_strands() const{ return max_strands_; }// how many strands to consider (if I see more than this number of strands, I ignore them
	void max_strands( core::Size const m ){ max_strands_ = m; }

	bool parallel() const{ return parallel_; }
	void parallel( bool const p ){ parallel_ = p; }

	void chain( core::Size const c ){ chain_ = c;}
	core::Size chain() const{ return chain_ ; }

	void sequence_separation( core::Size const c ){ sequence_separation_ = c;}
	core::Size sequence_separation() const{ return sequence_separation_; }

	void residues_to_align_on_pose( utility::vector1< core::Size > r ){ residues_to_align_on_pose_ = r; }
	void residues_to_align_on_template( utility::vector1< core::Size > r ){ residues_to_align_on_template_ = r; }
private:
	utility::vector1< core::Size > reference_positions( core::pose::Pose const & p ) const;
	core::Real distance_threshold_; // dflt 16;
	core::Size strand_length_, neighbors_, N_terminal_count_, max_strands_; //dflt 3, 6, 3, 10. Strand length: min number of contiguous residues in 'E' to count as a strand; neighbors: how many neighbors at distance_threshold does each N-terminal residue need to have to count in barrle; N-terminal count: How many residues should be aligned from each strand; max_strands: how many strands to align (specify 100000 to align all).
	bool odd_, even_; // dflt true, true; align even or odd (or both) strands
	core::pose::PoseOP template_pose_; //dflt NULL; pose to which to align
	core::Size stagger_; // dflt 0; see below
	bool parallel_; //dflt true; are the strands parallel or antiparallel?
	core::Size chain_; //dflt 1; on which chain to search? 0 means all chains
	core::Size sequence_separation_; // dflt 15aa; the minimal sequence separation between neighboring strands. Useful to eliminate short hairpins
	utility::vector1< core::Size > residues_to_align_on_pose_, residues_to_align_on_template_; //dflt empty; if specified short circuit computing which residues to align, instead aligning the residues specificed here directly.
};

/* explanation for stagger:
Stagger allows you to align a barrel to other barrels with a rotational stagger around the barrel principal axis. You can also align the barrel against itself.
For a beta barrel with m strands 1->m a stagger of n means that strands 1->(n) will be aligned to (m-n+1)->m and strands (n+1)->m will be aligned to strands 1->(m-n). E.g., For m=8 strands and stagger=5: (1->5):(4->8); (6->8):(1->3).
*/

} // simple_moves
} // protocols

#endif
