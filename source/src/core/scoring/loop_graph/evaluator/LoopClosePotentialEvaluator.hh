// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_LoopClosePotentialEvaluator_HH
#define INCLUDED_core_scoring_loop_graph_LoopClosePotentialEvaluator_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

	class LoopClosePotentialEvaluator: public utility::pointer::ReferenceCount {

	public:

		//constructor
		LoopClosePotentialEvaluator( core::Real const & loop_fixed_cost );

		//destructor
		~LoopClosePotentialEvaluator();

	public:

		core::Real const & loop_fixed_cost() const  { return loop_fixed_cost_; }

		void set_loop_closure_energy( core::Real const & setting ){ loop_closure_energy_ = setting; }
		core::Real loop_closure_energy() const { return loop_closure_energy_; }

		void set_involves_current_pose( bool const & setting ){ involves_current_pose_ = setting; }
		bool involves_current_pose() const { return involves_current_pose_; }

		void set_current_pose_takeoff_atom( core::id::AtomID const & setting ){ current_pose_takeoff_atom_ = setting; }
		core::id::AtomID current_pose_takeoff_atom() const { return current_pose_takeoff_atom_; }

		void set_current_pose_landing_atom( core::id::AtomID const & setting ){ current_pose_landing_atom_ = setting; }
		core::id::AtomID current_pose_landing_atom() const { return current_pose_landing_atom_; }

		void set_current_pose_takeoff_xyz( Vector const & setting ){ current_pose_takeoff_xyz_ = setting; }
		Vector current_pose_takeoff_xyz() const { return current_pose_takeoff_xyz_; }

		void set_current_pose_landing_xyz( Vector const & setting ){ current_pose_landing_xyz_ = setting; }
		Vector current_pose_landing_xyz() const { return current_pose_landing_xyz_; }

		virtual
		void
		get_f1_f2( Vector & f1, Vector & f2, bool const takeoff ) const = 0;


	private:

		core::Real loop_closure_energy_;
		bool involves_current_pose_;
		core::id::AtomID current_pose_takeoff_atom_;
		core::id::AtomID current_pose_landing_atom_;
		Vector current_pose_takeoff_xyz_;
		Vector current_pose_landing_xyz_;
		Real const loop_fixed_cost_; // offset to energy, defined by log V, were V is 'typical' restriction volume in Angstrom^3

	};

} //evaluator
} //loop_graph
} //scoring
} //core

#endif
