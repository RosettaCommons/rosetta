// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/MSDMover.hh
/// @brief Multistate design mover used for restrained multistate design
/// @brief Takes in multiple poses, applies residue linking constraints based on
/// @brief sequence of all input poses and runs a design submover that has been specified in the tag
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_MSDMover_hh
#define INCLUDED_protocols_simple_moves_MSDMover_hh

// Unit headers
#include <protocols/simple_moves/MSDMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/VectorPoseMover.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class MSDMover : public moves::VectorPoseMover {

public:
	/// @brief
	///  empty constructor fills values with the values
	///  read in from the commandline
	MSDMover();

	MSDMover( protocols::moves::MoverOP mover,
		utility::vector1< std::string > resfiles);

	~MSDMover() override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & ) override;


	void apply( Pose & pose ) override;

	void setup_mover ( Pose & pose );

	std::string get_name() const override;

	utility::vector1< core::scoring::constraints::ConstraintCOP >
	apply_linked_constraints( core::pose::Pose & pose );

	void parse_resfiles();

	utility::vector1< core::Size >
	parse_resfile ( core::pack::task::PackerTaskCOP design_task );

	moves::MoverOP design_mover();

	void design_mover( moves::MoverOP design_mover );

	utility::vector1< std::string > resfiles();

	void resfiles ( utility::vector1< std::string > resfiles );

	void weight( core::Real weight );

	core::Real weight();

	void update_packer_task ();

	void set_current_pose ( core::Size current_pose );

	utility::vector1< utility::vector1< core::Size > > res_links ();

private:
	std::string
	resfile_at ( core::Size index );

	protocols::moves::MoverOP design_mover_;
	utility::vector1< std::string > resfiles_;
	core::Real weight_;
	core::Size current_pose_;
	utility::vector1< utility::vector1< core::Size > > res_links_;

	bool debug_;
};


} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_MSDMover_HH
