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
#include <core/scoring/func/Func.fwd.hh>

namespace protocols {
namespace simple_moves {

class AddConstraintsToCurrentConformationMover : public moves::Mover {

public:
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;

	AddConstraintsToCurrentConformationMover();
	virtual ~AddConstraintsToCurrentConformationMover();

	virtual void apply( core::pose::Pose & );
	virtual std::string get_name() const;

	bool residue_to_constrain(Size const & i) const;

	void task_factory( TaskFactoryCOP tf );

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );

	/// @brief parse "task_operations" XML option (can be employed virtually by derived movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &,
		Filters_map const &,
		moves::Movers_map const &,
		Pose const & );


	bool       & use_distance_cst() { return use_distance_cst_; }
	bool       & CA_only() { return CA_only_; }
	bool       & bb_only() { return bb_only_; }
	bool       & inter_chain() { return inter_chain_; }
	core::Real & cst_weight() { return cst_weight_; }
	core::Real & max_distance() { return max_distance_; }
	core::Real & coord_dev() { return coord_dev_; }
	core::Real & bound_width() { return bound_width_; }
	core::Size & min_seq_sep() { return min_seq_sep_; }
	bool       const & use_distance_cst() const { return use_distance_cst_; }
	bool       const & CA_only() const { return CA_only_; }
	bool       const & bb_only() const { return bb_only_; }
	bool       const & inter_chain() const { return inter_chain_; }
	core::Real const & cst_weight() const { return cst_weight_; }
	core::Real const & max_distance() const { return max_distance_; }
	core::Real const & coord_dev() const { return coord_dev_; }
	core::Real const & bound_width() const { return bound_width_; }
	core::Size const & min_seq_sep() const { return min_seq_sep_; }


private:
	// pointers to data that are passed in
	core::pack::task::TaskFactoryCOP task_factory_;
	core::pack::task::PackerTaskOP task_;
	bool has_task_factory_;

	bool use_distance_cst_;
	bool CA_only_;
	bool bb_only_;
	bool inter_chain_;
	core::Real cst_weight_;
	core::Real max_distance_;
	core::Real coord_dev_;
	core::Real bound_width_;
	core::Size min_seq_sep_;

	core::scoring::func::FuncOP cc_func_;

};

} // moves
} // protocols

#endif
