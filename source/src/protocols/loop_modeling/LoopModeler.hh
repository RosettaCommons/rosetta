// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopModeler_HH
#define INCLUDED_protocols_loop_modeling_LoopModeler_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopModeler.fwd.hh>
#include <protocols/loop_modeling/LoopBuilder.fwd.hh>
#include <protocols/loop_modeling/LoopProtocol.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocols headers
#include <protocols/filters/Filter.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace loop_modeling {

class LoopModeler : public LoopMover {

public:

	/// @brief Default constructor.
	LoopModeler();

	/// @brief Default destructor.
	~LoopModeler();

	/// @brief Return the class name of this mover.
	string get_name() const { return "LoopModeler"; }

protected:

	bool do_apply(Pose & pose);

private:

	void convert_to_centroid(Pose & pose);
	void convert_to_fullatom(Pose & pose, Pose & original_pose);

public:

	void parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & data,
			protocols::filters::Filters_map const & filters,
			protocols::moves::Movers_map const & movers,
			Pose const & pose);

public:

	void enable_build_stage();
	void disable_build_stage();

	void enable_centroid_stage();
	void disable_centroid_stage();

	void enable_fullatom_stage();
	void disable_fullatom_stage();

	LoopBuilderOP build_stage();
	LoopProtocolOP centroid_stage();
	LoopProtocolOP fullatom_stage();

	void repack_everything_before_fullatom(bool setting);

	void setup_empty_config();
	void setup_kic_config();
	void setup_next_gen_kic_config();

	void add_shared_mover(LoopMoverOP mover);
	void add_shared_refiner(LoopMoverOP refiner);
	void add_shared_filter(protocols::filters::FilterOP filter);

	void clear_movers();
	void clear_refiners();
	void clear_movers_and_refiners();
	void mark_as_default();

private:

	LoopBuilderOP build_stage_;
	LoopProtocolOP centroid_stage_;
	LoopProtocolOP fullatom_stage_;

	bool is_build_stage_enabled_;
	bool is_centroid_stage_enabled_;
	bool is_fullatom_stage_enabled_;
	bool repack_everything_before_fullatom_;
};

}
}

#endif
