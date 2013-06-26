// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_simple_moves_AddConstraintsToCurrentConformationMover_hh
#define INCLUDED_protocols_simple_moves_AddConstraintsToCurrentConformationMover_hh

#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>


namespace protocols {
namespace simple_moves {

class AddConstraintsToCurrentConformationMover : public moves::Mover {

public:
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;

	AddConstraintsToCurrentConformationMover();
	virtual ~AddConstraintsToCurrentConformationMover();

	virtual void apply( core::pose::Pose & );
	virtual std::string get_name() const;

	void task_factory( TaskFactoryCOP tf );

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagPtr const, moves::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );

	///@brief parse "task_operations" XML option (can be employed virtually by derived movers)
	virtual void parse_task_operations(
									   TagPtr const,
									   moves::DataMap const &,
									   Filters_map const &,
									   moves::Movers_map const &,
									   Pose const & );

private:
	// pointers to data that are passed in
	core::pack::task::TaskFactoryCOP task_factory_;
  core::pack::task::PackerTaskOP task_;
	bool has_task_factory_;

	bool use_distance_cst_;
	bool CA_only_;
	core::Real cst_weight_;
	core::Real max_distance_;
	core::Real coord_dev_;
	core::Real bound_width_;
	core::Size min_seq_sep_;

};

} // moves
} // protocols

#endif
