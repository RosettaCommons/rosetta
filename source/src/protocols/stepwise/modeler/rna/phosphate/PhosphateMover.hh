// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/phosphate/PhosphateMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_rna_phosphate_PhosphateMover_HH
#define INCLUDED_protocols_stepwise_modeler_rna_phosphate_PhosphateMover_HH

#include <protocols/stepwise/modeler/rna/phosphate/PhosphateMover.fwd.hh>
#include <protocols/stepwise/modeler/rna/phosphate/PhosphateMove.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace phosphate {


	///////////////////////////////////////////////////////////////
	class PhosphateMover: public protocols::moves::Mover {

	public:

		//constructor
		PhosphateMover( Size const sample_res,
										PhosphateTerminus which_terminus_,
										scoring::ScoreFunctionCOP scorefxn );

		//constructor
		PhosphateMover( PhosphateMove const phosphate_move,
										scoring::ScoreFunctionCOP scorefxn );

		//destructor
		~PhosphateMover();

		virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

	public:

		void
		screen_phosphate( pose::Pose & pose );

		bool instantiated_phosphate() const{ return instantiated_phosphate_; }

		void set_force_phosphate_instantiation( bool const & setting ){ force_phosphate_instantiation_ = setting; }
		bool force_phosphate_instantiation() const { return force_phosphate_instantiation_; }

	private:

		void
		initialize_variables();

		void
		setup_variants_and_free_pose_for_terminal_phosphate( pose::Pose & pose  );

		void
		setup_variants_and_free_pose_for_five_prime_phosphate( pose::Pose & pose );

		void
		setup_variants_and_free_pose_for_three_prime_phosphate( pose::Pose & pose );

		void
		screen_five_prime_phosphate( pose::Pose & pose );

		void
		screen_three_prime_phosphate( pose::Pose & pose );

		void
		setup_atom_and_neighbor_list( pose::Pose & pose );

		bool
		check_phosphate_contacts_donor( pose::Pose & pose ) const;

		bool
		pass_clash_check( std::string const atom_name,
											Size const n,
											pose::Pose & pose );

		void
		apply_Aform_torsions_to_five_prime_phosphate( pose::Pose & pose, Size const sample_res ) const;

		void
		apply_Aform_torsions_to_three_prime_phosphate( pose::Pose & pose, Size const sample_res ) const;

	private:

		PhosphateMove const phosphate_move_;
		scoring::ScoreFunctionCOP scorefxn_;

		bool do_screening_;
		bool screen_for_donor_contact_;
		bool instantiated_phosphate_;
		bool force_phosphate_instantiation_;

		pose::PoseOP pose_free_;

		utility::vector1< Vector > donor_atom_xyz_list_;
		utility::vector1< Vector > donor_base_atom_xyz_list_;
		utility::vector1< Size > neighbor_copy_dofs_;
		Size op1_atom_idx_, op2_atom_idx_;
		Size number_score_calls_;

		core::chemical::rna::RNA_FittedTorsionInfo const torsion_info_;

	};

} //phosphate
} //rna
} //modeler
} //stepwise
} //protocols

#endif
