// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/MergePDBMover.hh
/// @brief This class will allign & combine parts of the pdb.
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_MergePDBMover_hh
#define INCLUDED_protocols_simple_moves_MergePDBMover_hh


#include <protocols/simple_moves/MergePDBMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace simple_moves {


class MergePDBMover : public moves::Mover {
public:
	MergePDBMover();

	void determine_overlap(Pose const pose, Size & start_overlap_parent, Size & end_overlap_parent, Size & start_overlap_pose, Size & end_overlap_pose);
	void generate_overlap(Pose & pose, Size start_overlap_parent, Size end_overlap_parent, Size start_overlap_pose, Size end_overlap_pose);
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

private:
	core::pose::PoseOP parent_pose_;
	std::string overlap_location_pose_;
	Size overlap_range_;
	Size overlap_max_length_;
};

} // simple_moves
} // protocols

#endif

