// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/rna/RNA_ResampleMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_RNA_ResampleMover_HH
#define INCLUDED_protocols_stepwise_monte_carlo_RNA_ResampleMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_ResampleMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/SWA_Move.fwd.hh>
#include <protocols/stepwise/monte_carlo/SWA_MoveSelector.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Modeler.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	class RNA_ResampleMover: public protocols::moves::Mover {

	public:

		//constructor
		RNA_ResampleMover( protocols::stepwise::enumerate::rna::StepWiseRNA_ModelerOP stepwise_rna_modeler );

		//destructor
		~RNA_ResampleMover();

	public:

		/// @brief Apply the minimizer to one pose
		virtual void apply( core::pose::Pose & pose_to_visualize );
		virtual std::string get_name() const;

		bool
		apply( core::pose::Pose & pose,
					 std::string & move_type );

		bool
		apply( core::pose::Pose & pose,
					 SWA_Move & swa_move );

		bool
		apply( core::pose::Pose & pose,
					 SWA_Move const & swa_move,
					 std::string & move_type );


		void set_just_min_after_mutation_frequency( core::Real const & setting ){ just_min_after_mutation_frequency_ = setting; }
		core::Real just_min_after_mutation_frequency() const{ return just_min_after_mutation_frequency_; }

		void set_allow_internal_moves( bool const & setting ){ allow_internal_moves_ = setting; }
		bool allow_internal_moves() const{ return allow_internal_moves_; }

		void set_minimize_single_res( bool const & setting ){ minimize_single_res_ = setting; }
		bool minimize_single_res() const{ return minimize_single_res_; }

	private:

		protocols::stepwise::enumerate::rna::StepWiseRNA_ModelerOP stepwise_rna_modeler_;
		protocols::stepwise::monte_carlo::SWA_MoveSelectorOP swa_move_selector_;

		core::Real just_min_after_mutation_frequency_;
		bool allow_internal_moves_;
		bool minimize_single_res_;
		core::Real constraint_x0_;
		core::Real constraint_tol_;

	};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
