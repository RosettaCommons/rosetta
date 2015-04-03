// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/loop_graph/LoopGraph.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_LoopGraph_HH
#define INCLUDED_core_scoring_loop_graph_LoopGraph_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/loop_graph/LoopGraph.fwd.hh>
#include <core/scoring/loop_graph/Loop.fwd.hh>
#ifdef WIN32
#include <core/scoring/loop_graph/Loop.hh>
#include <core/scoring/loop_graph/LoopCycle.hh>
#endif
#include <core/scoring/loop_graph/LoopCycle.fwd.hh>
#include <core/scoring/loop_graph/LoopScoreInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>

#include <map>

namespace core {
namespace scoring {
namespace loop_graph {

	class LoopGraph: public utility::pointer::ReferenceCount {

	public:

	//constructor
	LoopGraph();

	//destructor
	~LoopGraph();

	public:

		void update( pose::Pose & pose, bool const verbose = false );

		LoopScoreInfoOP loop_score_info( Size const n ) const;

		Size num_loops() const{ return current_pose_loop_score_info_.size(); }

		Real total_energy() const{ return total_energy_; }

		// helper functions for computing loop properties, missing residues, etc.
		void update_loops( core::pose::Pose const & pose );

		void
		update_loops_and_cycles( utility::vector1< Size > const & pose_domain_map,
														 utility::vector1< Size > const & cutpoint_open );

		Size nmissing( pose::Pose const & pose ) const;

		utility::vector1< char >
		missing_residues( pose::Pose const & pose ) const;

		//////////////////////////////////////////////////////////////////
		/// @details Returns a vector of loop_suites
		/// @brief   If include_free_loops is turned off, loops with 
		/// takeoff_domain == 0 or landing_domain == 0 will not be included
		utility::vector1< utility::vector1< Size > >
		loop_suites( bool include_free_loops = true ) const;

		void set_error_out_on_complex_cycles( bool const & setting ){ error_out_on_complex_cycles_ = setting; }
		bool has_just_simple_cycles() const { return has_just_simple_cycles_; }

		void
		check_loop_cycles_are_disjoint( bool const verbose = false );

	private:

		void
		get_loop_atom( Size const & res,
									 core::pose::Pose const & pose,
									 bool const takeoff /* as opposed to landing */,
									 id::AtomID & atom_id,
									 Vector & xyz );

		void
		figure_out_loop_cycles();

		void
		look_for_cycles_recursively( Size const current_domain,
																 utility::vector1< Size > parent_domains_in,
																 utility::vector1< Loop > loops_so_far_in );

		bool
		check_disjoint( LoopCycle loop_cycle1, LoopCycle loop_cycle2 ) const;

		void
		update_loops( utility::vector1< Size > const & pose_domain_map,
														 utility::vector1< Size > const & cutpoint_open );

		void
		check_for_unexpected_cutpoints( pose::Pose const & pose ) const;

	private:

		utility::vector1< Loop > loops_;
		utility::vector1< LoopCycle > loop_cycles_;
		utility::vector1< LoopScoreInfoOP > current_pose_loop_score_info_;

		std::map< Size, utility::vector1< Size > > loops_from_domain_;
		utility::vector1< bool > loop_visited_;
		utility::vector1< bool > domain_visited_;

		Real const rna_gaussian_variance_per_residue_;
		Real const protein_gaussian_variance_per_residue_;
		Real const loop_fixed_cost_;

		Real total_energy_;
		bool error_out_on_complex_cycles_;
		bool has_just_simple_cycles_;

	};

} //loop_graph
} //scoring
} //core

#endif
