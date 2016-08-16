// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/relax/AtomCoordinateCstMover.hh
/// @brief Coordinate constraints for relax, etc.
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_relax_AtomCoordinateCstMover_hh
#define INCLUDED_protocols_relax_AtomCoordinateCstMover_hh

// Unit headers
#include <protocols/relax/AtomCoordinateCstMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <protocols/loops/Loops.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/types.hh>

// C++ headers
#include <string>

namespace protocols {
namespace relax {

class AtomCoordinateCstMover : public moves::Mover {
public:
	typedef moves::Mover parent;

public:

	AtomCoordinateCstMover();
	~AtomCoordinateCstMover();

	virtual void
	apply( core::pose::Pose & pose );

	std::string get_name() const { return "AtomCoordinateCstMover"; }

	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual protocols::moves::MoverOP clone() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	virtual void set_refstruct( core::pose::PoseCOP ref );

	void cst_sd( core::Real sd ) { cst_sd_ = sd; }
	core::Real cst_sd() const { return cst_sd_; }

	void bounded( bool bound ) { bounded_ = bound; }
	bool bounded( ) const { return bounded_; }

	void cst_width( core::Real width ) { cst_width_ = width; }
	core::Real cst_width() const { return cst_width_; }

	void cst_sidechain( bool sidechains ) { cst_sidechain_ = sidechains; }
	bool cst_sidechain( ) const { return cst_sidechain_; }

	void ambiguous_hnq( bool flip_hnq ) { amb_hnq_ = flip_hnq; }
	bool ambiguous_hnq( ) const { return amb_hnq_; }

	void set_loop_segments( protocols::loops::LoopsCOP loops);
	void set_task_segments( core::pack::task::TaskFactoryCOP taskfactory);

private:
	core::select::residue_selector::ResidueSubset
	compute_residue_subset( core::pose::Pose const & pose ) const;

	core::pose::PoseOP
	get_constraint_target_pose( core::pose::Pose const & pose ) const;

	core::id::SequenceMapping
	generate_seqmap( core::pose::Pose const & pose, core::pose::Pose const & constraint_target_pose ) const;

	core::scoring::constraints::ConstraintCOPs
	generate_constraints( core::pose::Pose const & pose );

private:
	/// @brief If set, the pose to make the constraints to
	core::pose::PoseCOP refpose_;

	/// @brief The strenght (standard deviation) of the constraints
	core::Real cst_sd_;

	/// @brief Use bounded constraints instead of harmonic
	bool bounded_;

	/// @brief The width for bounded constraints
	core::Real cst_width_;

	/// @brief Also constrain sidechains?
	bool cst_sidechain_;

	/// @brief When making sidechain constraints, should we allow for flipping HNQ residue sidechains?
	bool amb_hnq_;

	/// @brief A loop definition of constrained segments
	protocols::loops::LoopsCOP loop_segments_;

	/// @brief A task definition of constrained segments
	core::pack::task::TaskFactoryCOP task_segments_;
};

}
} // protocols

#endif
