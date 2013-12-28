// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/rna/TransientCutpointHandler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_TransientCutpointHandler_HH
#define INCLUDED_protocols_stepwise_TransientCutpointHandler_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/monte_carlo/rna/TransientCutpointHandler.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	class TransientCutpointHandler: public utility::pointer::ReferenceCount {

	public:

	//constructor
		TransientCutpointHandler( core::Size const sample_res );

		TransientCutpointHandler( core::Size const sample_suite,
															core::Size const cutpoint_suite );

	//destructor
	~TransientCutpointHandler();

	public:

		void put_in_cutpoints( core::pose::Pose & pose );

		void take_out_cutpoints( core::pose::Pose & pose );

		void set_minimize_res( utility::vector1< Size > const & setting ){ minimize_res_ = setting; }

		void set_move_jump_points_away( bool const & setting ){ move_jump_points_away_ = setting; }
		bool move_jump_points_away() const{ return move_jump_points_away_; }


	private:

		void prepare_fold_tree_for_erraser( core::pose::Pose & pose );

		core::Size const sample_suite_, cutpoint_suite_;
		bool move_jump_points_away_;

		utility::vector1< core::Size > fixed_res_;
		utility::vector1< core::Size > minimize_res_;
		core::kinematics::FoldTree fold_tree_save_;

	};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
