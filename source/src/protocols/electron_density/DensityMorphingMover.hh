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

#ifndef INCLUDED_protocols_electron_density_DensityMorphingMover_hh
#define INCLUDED_protocols_electron_density_DensityMorphingMover_hh

#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/id/AtomID.hh>

namespace protocols {
namespace electron_density {

class DensityMorphingMover : public moves::Mover {

public:

	DensityMorphingMover();
	virtual ~DensityMorphingMover();

	void init();

	utility::vector1<core::id::AtomID> collect_fragment_atom_ids(core::pose::Pose const & pose, core::Size ires, core::Size extend_residues);

	virtual void apply( core::pose::Pose & );

	//void task_factory( core::pack::task::TaskFactoryCOP tf );

	std::string get_name() const { return "DensityMorphingMover"; }

	moves::MoverOP clone() const { return new DensityMorphingMover( *this ); }
	moves::MoverOP fresh_instance() const { return new DensityMorphingMover; }

	virtual void
	parse_my_tag( TagCOP const, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );

	///@brief parse "task_operations" XML option (can be employed virtually by derived movers)
	/*
	virtual void parse_task_operations(
									   TagCOP const,
									   basic::datacache::DataMap const &,
									   Filters_map const &,
									   moves::Movers_map const &,
									   Pose const & );
	*/
	
private:
	core::Size nres_;
	core::Size frag_size_;
	core::Real search_radius_;
	core::Real cst_weight_;
	core::Real bound_width_;
	core::Real coord_dev_factor_;

};

} // moves
} // protocols

#endif
