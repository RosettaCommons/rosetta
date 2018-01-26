// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyLoophashAssemblyMover.hh
///
/// @brief A Mover that uses loophash to find fragments that can
/// bridge a gap with minimal modifications to the original pose.
/// @author Tim Jacobs

#ifdef NOT_IN_SCONS_DEPRECATED

#ifndef INCLUDED_protocols_legacy_sewing_sampling_LegacyLoophashAssemblyMover_HH
#define INCLUDED_protocols_legacy_sewing_sampling_LegacyLoophashAssemblyMover_HH

//Unit headers
#include <protocols/legacy_sewing/sampling/LegacyLoophashAssemblyMover.fwd.hh>

//Package headers
#include <protocols/legacy_sewing/sampling/LegacyAssemblyMover.hh>
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Protocol headers
#include <core/pose/Pose.hh>

#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>

#include <protocols/loops/Loops.hh>


namespace protocols {
namespace legacy_sewing  {

class LegacyLoophashAssemblyMover : public LegacyAssemblyMover {

public:

	typedef utility::vector1< std::pair< protocols::loophash::BackboneSegment, std::string > > BackboneSegments;

	LegacyLoophashAssemblyMover();

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	std::string
	get_name() const;

	//void init();

	///@brief override complete assembly to include building
	///of Loophash segments
	virtual
	bool
	complete_assembly(
		AssemblyOP & assembly
	);

	///@brief override statistics output
	///to include some loop-related output
	virtual
	void
	output_stats(
		AssemblyOP const & assembly,
		core::pose::Pose & pose
	);


	bool
	rearrange_assembly(
		AssemblyOP & assembly
	) const;

	///@brief count the number of loophash segments between unconnected jumps in the assembly
	core::Size
	count_loophash_fragments(
		AssemblyOP const assembly,
		core::pose::Pose const & pose
	) const;

	///@brief try to add a loophash segments to the assembly
	protocols::loops::Loops
	add_loophash_segments(
		AssemblyOP & assembly,
		core::pose::Pose & pose
	) const;

	///@brief try to add a signle loophash segment to the pose at the anchor residue
	protocols::loops::Loop
	add_single_loop(
		core::pose::Pose & pose,
		core::Size loop_anchor,
		core::Size n_segment_start,
		core::Size c_segment_end
	) const;

	///@brief Get the backbone segments between loop_anchor and loop_anchor+1
	BackboneSegments
	get_backbone_segments(
		core::pose::Pose & pose,
		core::Size loophash_fragment_start,
		core::Size loophash_fragment_end
	) const;

	///@brief Build residues and cart-min close the given loophash backbone segment
	std::pair<core::pose::Pose, core::Real>
	build_loop_pose(
		core::pose::Pose const & pose,
		core::Size loop_anchor,
		core::Size n_segment_start,
		core::Size c_segment_end,
		core::Size loophash_fragment_start,
		core::Size loophash_fragment_end,
		protocols::loophash::BackboneSegment bb_seg,
		std::string loop_sequence
	) const;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose
	);

	///@brief Remove all backbone segments not within the torsion RMS
	//limits for the flanking regions of the loop
	void
	trim_bb_segs(
		core::pose::Pose const & pose,
		core::Size loop_anchor,
		core::Size n_overlapping_n,
		core::Size n_overlapping_c,
		BackboneSegments & bb_segs
	) const;

private:

	protocols::loophash::LoopHashLibraryCOP lh_library_;
	protocols::loops::Loops built_loops_;
	core::Size max_loop_segments_to_build_;
	core::Real max_loop_distance_;

	core::Size min_flanking_residues_;
	core::Size max_flanking_residues_;

	bool check_flanking_rms_;
	core::Real flanking_rms_cutoff_;

	core::Real max_insertion_rms_;

	core::scoring::ScoreFunctionOP loop_scorefxn_;
	core::scoring::ScoreFunctionOP loop_refine_scorefxn_;

};

} //legacy_sewing
} //protocols

#endif

#endif
