// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SimpleThreadingMover.cc
/// @brief Very Simple class for threading a regional sequence onto a structure
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_SimpleThreadingMover_hh
#define INCLUDED_protocols_simple_moves_SimpleThreadingMover_hh

#include <protocols/simple_moves/SimpleThreadingMover.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


// Forward
namespace protocols {
namespace simple_moves {

///@brief
/// This mover functions to thread the sequence of a region onto the given pose.  Nothing fancy here.
/// For more fancy things see protocols/comparative_modeling.
///
///@details
/// It does the threading by allowing the task to only enable these residues and then does a repacking. Optionally repack neighbors so we save one more step.
/// A sequence is just a string, additional '-' charactors denote to skip this position in the thread.
/// Default is 5 rounds of packing.
///
class SimpleThreadingMover : public protocols::moves::Mover {

public :

	SimpleThreadingMover();
	SimpleThreadingMover(std::string thread_sequence, core::Size start_position);

	SimpleThreadingMover(SimpleThreadingMover const & src);

	~SimpleThreadingMover() override;


	///@brief Set the sequence to thread onto the structure used in apply and where to start.
	///@details Can have '-' charactors in sequence to denote a gap in the threaded sequence.
	void
	set_sequence(std::string thread_sequence, core::Size start_position);

	///@brief Pack the neighbor residues?
	void
	set_pack_neighbors(bool pack_neighbors);

	///@brief Set the packing distance for neighbor pack.
	void
	set_neighbor_distance(core::Real neighbor_dis);

	///////////////

	bool
	get_pack_neighbors() const;

	core::Real
	get_neighbor_distance() const;

	///@brief Set the scorefunction used for packing.
	void
	set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn);

	///@brief Set the number of pack rounds.
	void
	set_pack_rounds(core::Size pack_rounds);

	void
	apply(core::pose::Pose & pose) override;


public:

	// XRW TEMP  std::string
	// XRW TEMP  get_name() const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose) override;

	protocols::moves::MoverOP
	clone() const override;

	//SimpleThreadingMover & operator=( SimpleThreadingMover const & src);

	moves::MoverOP fresh_instance() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	void
	set_defaults();

private:

	core::Size start_position_;
	std::string thread_sequence_;

	bool pack_neighbors_;
	core::Real neighbor_dis_;

	core::scoring::ScoreFunctionCOP scorefxn_;

	std::string parsed_position_; //enables pose-length changes after construction of mover.
	bool skip_unknown_mutant_;

	core::Size pack_rounds_;

};

}//simple_moves
}//protocols



#endif //INCLUDED_protocols_simple_moves_SimpleThreadingMover_hh







